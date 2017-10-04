/* contains c_bijrhsmarvex */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++) 
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Sqr(a) ((a)*(a))
#define print2d(name,var,im,jm) {printout("normal","%s %d %d ",name,im,jm); Loop(ix,0,im) Loop(jx,0,jm) printout("normal"," %lg",var[ix][jx]); printout("normal","\n"); }
#define print1d(name,var,im) {printout("normal","%s %d ",name,im); Loop(ix,0,im) printout("normal"," %lg",var[ix]); printout("normal","\n");}

/* pressure strain functions */
#define fg0(gp,gl) ((max(0,((gl)-(gp))/(gl)))*(max(0,((gl)-(gp))/(gl))))
#define fg1(gp,gl) ((max(0,((gp)-(gl))/(1.-(gl))))*(max(0,((gp)-(gl))/(1.-(gl)))))
#define flin(gp,gl,gh) (max(0,(gh)-(gl)-max(0,(gp)-(gl)))/((gh)-(gl)))
#define fg3(gp,gl,gh)  (1.+flin(gp,gl,gh)* flin(gp,gl,gh)*(-3.+2.* flin(gp,gl,gh)))

/*------------------  MARV  model-------------------*/
/*  with options of including (or omitting) 'rvs' */
/*  with optional added c5 and/or c7 term for non-linear instabilities */
/*  with optional pkdktot in place of local pkdk calc here */

void ec_bijrhsmarvex(FILE *fpin, FILE *fprint)
/* in vki notes
 bij eqn source term =
 -beqc b[ijcen] +beqm[ij] + beqd[ij][nm] db[nmcen]
 in this routine include  beqc in beqd and beqm (beqc is temporary in this routine)
 to give
 bij eqn source term = beqm[ij] + beqd[ij][nm] db[nmcen]
 this routine gives  beqm[noindfwpts*6], beqd[noindfwpts*6*6]
 also gives intermediate arrays (c1 ...) in beqex if bug>0 
 and between the points values for pkdk (Prod/k) and gbij (MARV g parameter)
 
 bijsourcemarv revised to give rhsbij full rhs including coef 0 contributions
 also set cambij[noindfwpts] = max of sum of + or - coefficients in beqd
 beqm same as rhsbij, where call beqm source only rhsbij with added coef 0 contribution 
 */

{ 
  int *i4d,*wherep,*noindfwpts,*wherefw,*whoisfw;   /* need */
  char *cltp;
  double *u[3],*rho,*q,*om,*vlam,*bij;
  int *cn,**ci; double **cc;
  int *itrange; double *zrotation; /* use if available */
  double *walldist=0; /* needed for rotating frame, walls are currently assumed to rotate */
  double *beqm, *cambij,*rhsbij,*beqd, *pkdk, *gbij;  /* set bij eq coef, turbulence P/k, g */
  
  double *beqex=0; /* optionally set mid-value array with extras, c1 etc. (see print) */
  double *beqc; /* temporary */
  int ibeqex=12, lpde=0, lrt=1, lgeff=2, lfw=3, lfb=4, lc1=5, lc3=6, lc4=7, lc6=8,lfr=9,lcapg=10, lc1homo=11;
  char Nbeqex[]=" pde~,rt~,geff,fw,fb,c1,c3,c4,c6,fr,capg,c1homo";
  int bug=1;   /* *******   set bug ****** */
  
  int i,j,k,n,m,ip[4],i4dp[4],i4dm[4],its,ite,mc[8],mmid,ia[3],ii[4],iall,iallm,ieq,nf,iter;
  int ix,jx;
  int two[3]={2,2,2};
  double cave[8],volc,vol[8],grad[8][3];
  double bm[3][3],dudx[3][3],sij[3][3],wij[3][3],wijr[3][3];
  double rot[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double bb,bbb,gradg[3],g,gradq[3],omm,qm,rhom,ss,bmnsmn,bmnbmksnk;
  double pdk,rt,epst,eps;
  double geff,capg,c1inf,fc1,c1;
  double c3fdfz,c3fdfs,c3rdtz,c3rdts,c3fdf,c3rdt,c3;
  double c4fdfz,c4fdfs,c4rdtz,c4rdts,c4fdf,c4rdt,c4;
  double c6fdfz,c6fdfs,c6rdtz,c6rdts,c6fdf,c6rdt,c6;
  double fw,capb,fb,denom,bparm0,bparm1,bpold,bpnew,bplo,bphi, frot;
  double c1xp,c1xm,c3x,c4x;
  double tbeqc,tbeqm[6],tbeqd[6][6];
  double dp,dm,dd;
  
  int ij2n[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
  int n2ii[6]    = {0,1,2,0,0,1};
  int n2jj[6]    = {0,1,2,1,2,2};
  
  double *c5=0, *c7=0, *pkdktot=0;  /* for additional options */
  char *namemodel,*namec5,*namec7,*namepkdktot; /* input parameters */
  char modr='n',modv='n';
  double acapg;

  namemodel=readname(fpin);
  namec5=readname(fpin);
  namec7=readname(fpin);
  namepkdktot=readname(fpin); 
  printout("normal","using %s model + preset c5=%s,  c7=%s,  pkdktot=%s\n",
          namemodel,namec5,namec7,namepkdktot);

  if (namec5[0]!='\0') c5=(double *)need(namec5); 
  if (namec7[0]!='\0') c7=(double *)need(namec7);  
  if (namepkdktot[0]!='\0') pkdktot=(double *)need(namepkdktot);
  acapg=0.3;
  for (i=0;i<10;i++)
  {
    if (namemodel[i]=='\0') break;
    if (namemodel[i]=='r' || namemodel[i]=='R') modr='y';
    if (namemodel[i]=='v' || namemodel[i]=='V') modv='y';
    if (namemodel[i]=='s' || namemodel[i]=='S') acapg=0.244;
  }
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallm=Prod4(i4dm);
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  noindfwpts=(int *)need("noindfwpts");
  wherefw=(int *)need("wherefw");
  whoisfw=(int *)need("whoisfw");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  q=(double *)need("qturb");
  om=(double *)need("omturb");
  bij=(double *)need("bij");
  vlam=(double *)need("vlam");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1; }
  zrotation=(double *)find("zrotation");
  if (zrotation>0) 
  {  
    rot[0][1]=-zrotation[0]; 
    rot[1][0]=zrotation[0];
    walldist=(double *)need("walldist");
  }
  geomcinit();
  geom8init();
  
  beqm=(double*)createarray("rhsbij",noindfwpts[0]*6,'d',0);
  rhsbij=beqm;
  cambij=(double*)createarray("cambij",noindfwpts[0],'d',0);
  beqd=(double*)createarray("beqd",noindfwpts[0]*6*6,'d',0);
  pkdk=(double*)createarray("pkdk",iallm,'d',0);
  gbij=(double*)createarray("gbij",iallm,'d',0);
  Loop(i,iallm/i4d[3]*its,iallm/i4d[3]*(ite+1)) { pkdk[i]=0; gbij[i]=1;}
  
  beqc=(double *)tmalloca(noindfwpts[0],'d');
  if (bug>0)
  { 
    beqex=(double *)find("beqex");
    if (beqex==0) beqex=(double*)createarray("beqex",iallm*ibeqex,'d',0);
    Loop(i,0,iallm*ibeqex) beqex[i]=0;
    printout("normal"," saving %s in beqex\n",Nbeqex);
  }
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])
  {  
    beqc[i]=0; 
    Loop(j,0,6) beqm[i+j*noindfwpts[0]]=0;
    Loop(j,0,36) beqd[i+j*noindfwpts[0]]=0; 
  }
  /* loop over contiuity control volumes */
  Loop(ip[3],its,ite+1) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  if (wherep[In4(ip,i4dp)]>=0)   /* do only for non-zero cont c.v. */
  { 
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));  /* viscosity index */
    if (bug>1) printout("normal","ip %d %d %d %d mmid %d\n",ip[0],ip[1],ip[2],ip[3],mmid); 
    volc=geomcgradvol(ip,grad);   /* geometry arrays */
    if (bug>2) print2d("grad",grad,8,3);
    geomcvave(cave,ip);
    if (bug>2) print1d("cave",cave,8);
    geom8vol(vol,ip);
    if (bug>2) print1d("vol",vol,8);
    
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }         
    Loop(i,0,3) Loop(j,0,3)   /*  velocity gradient */
    { 
      dudx[i][j]=0;
      Loop(k,0,8) dudx[i][j] += grad[k][j]*u[i][mc[k]];
    }
    fixdudx(dudx);    /* fix dudx for incompressible continuity */
    if (bug>2) print2d("dudx",dudx,3,3);
    Loop(i,0,3) Loop(k,0,3) /* strain rate and vorticity */
    { 
      sij[i][k]=.5*(dudx[i][k]+dudx[k][i]);
      wijr[i][k]=.5*(dudx[i][k]-dudx[k][i]);
      wij[i][k]=wijr[i][k]+rot[i][k];
    }
    if (bug>2) {print2d("sij",sij,3,3); print2d("wij",wij,3,3);}
    ss=0;  /* strain rate */
    Loop(i,0,3) Loop(k,0,3) ss += sij[i][k]*sij[i][k];
    ss*= 2;
    if (bug>2) printout("normal","ss %lg\n",ss);
    Loop(j,0,3)   /* grad q */
    { 
      gradq[j]=0;
      Loop(k,0,8) gradq[j] +=grad[k][j]*q[mc[k]];
    }
    if (bug>2) print1d("gradq",gradq,3);
    Loop(i,0,3) gradg[i]=0;	  /* grad g using on the points g*/
    Loop(k,0,8)     
    { 
      bb=0; bbb=0;
      Loop(i,0,3) Loop(j,0,3)
      { 
        bb+=bij[mc[k]+iall*ij2n[i][j]]*bij[mc[k]+iall*ij2n[i][j]];
        Loop(n,0,3) bbb+=bij[mc[k]+iall*ij2n[i][j]]
        *bij[mc[k]+iall*ij2n[j][n]]*bij[mc[k]+iall*ij2n[n][i]];
      }
      g=1.-4.5*bb+9.*bbb;
      Loop(i,0,3) gradg[i]+=g*grad[k][i];
    }
    if (bug>2) print1d("gradg",gradg,3);
    Loop(i,0,3) Loop(j,0,3) bm[i][j]=0;  /* average values rho q om bij */
    omm=0; qm=0; rhom=0;
    Loop(k,0,8)  
    { 
      rhom+=cave[k]*rho[mc[k]];
      /*   omm+=cave[k]/om[mc[k]]; */ /* assume 1/om varies linearly */
      omm+=cave[k]*log(om[mc[k]]);
      qm+=cave[k]*q[mc[k]];
      Loop(i,0,3) Loop(j,0,3) bm[i][j] += cave[k]*bij[mc[k]+iall*ij2n[i][j]];
    }
    /* omm=1./omm; */
    omm=exp(omm);     
    if (bug>2) printout("normal","rho om q %lg %lg %lg\n",rhom,omm,qm);
    if (bug>2) print2d("bm",bm,3,3);
    bb=0; bbb=0;  /* redo g,bb,bbb using mid bij's */
    Loop(i,0,3) Loop(j,0,3)
    { 
      bb+=bm[i][j]*bm[i][j];
      Loop(n,0,3) bbb+=bm[i][j]*bm[j][n]*bm[n][i];
    }
    g=1.-4.5*bb+9.*bbb;
    gbij[mmid]=g;
    if (g<0) g=0; else if (g>1) g=1;
    if (bug>2) printout("normal","g %lg\n",g);
    bmnsmn=0.;  bmnbmksnk=0.;  /* b-s products*/
    Loop(n,0,3) Loop(i,0,3)
    { 
      bmnsmn += bm[i][n] * sij[i][n];
      Loop(k,0,3) bmnbmksnk += bm[i][n] * bm[i][k]* sij[k][n];
    }
    pdk=-2.*bmnsmn;    /* P/k Rt dissipation */
    pkdk[mmid]=pdk;   /* save for use elsewhere */
    if (pkdktot > 0) pdk=pkdktot[mmid];  /* use alt if exists */
    rt=0; epst=0;
    if (omm>0) 
    { 
      rt=rhom*qm*qm/(omm*vlam[mmid]);
      epst=omm*qm*qm;
    }
    eps=epst+2*(vlam[mmid]/rhom)*(gradq[0]*gradq[0]+gradq[1]*gradq[1]+gradq[2]*gradq[2]);
    if (bug>0) {beqex[mmid+iallm*lpde]=pdk/omm; beqex[mmid+iallm*lrt]=rt; }
    
    /* MARV model parameters */
    fw=0;  /*  fw parameter */
    if (ss>0)
    { 
      Loop(i,0,3) Loop(j,0,3) Loop(k,0,3) Loop (n,0,3)
      fw+=(wij[i][j]*sij[j][k]/ss)*(wij[k][n]*sij[n][i]/ss);
      fw *=8;
      if (fw<0) fw=0; else if (fw>1) fw=1;
      fw=sqrt(fw);
    }
    /* no V option */
    if (modv=='n') fw=1;
    
    geff=g;    /* effective g for MARV  model and c1 */
    capg=0; 
    if (g>.0001)
    { 	
      if (qm>0 && omm>0) capg=(qm/omm)*(qm/omm)
        *(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);
      geff=g*g/(g+acapg*capg*sqrt(acapg*capg));
    }
    if (bug>0) 
    {
      beqex[mmid+iallm*lgeff]=geff; 
      beqex[mmid+iallm*lfw]=fw;
      beqex[mmid+iallm*lcapg]=capg;
    }
    fc1=.000003*rt*rt+.045*sqrt(rt);	
    if (fc1>.001) fc1=1.-exp(-fc1);
    /* save c1 without geff for interest */
    if (bug>0)
    { 
      c1inf=8.+(2.-8.)*fg0(g,.5)+(2.5-8.)*fg1(g,.7); 
      c1=2.+(c1inf-2.)*fc1;
      beqex[mmid+iallm*lc1homo]=c1;
    }
    c1inf=8.+(2.-8.)*fg0(geff,.5)+(2.5-8.)*fg1(geff,.7); 
    c1=2.+(c1inf-2.)*fc1;
    if (bug>0) beqex[mmid+iallm*lc1]=c1;
    /* c3 c4 c6 for MARV model */
    c3fdfz=.44 -.44*Sqr(max(0,1.-sqrt(4.*g)))+.36*fg3(g, .7, 1.);
    c3fdfs=.15*fg3(g, .05, .5)+.65*fg3(g, .7, 1.);
    c3rdtz=.6+.2*g;
    c3rdts=.1-.1*Sqr(max(0, 1-sqrt(g/.15)))+.7*g;
    
    c4fdfz=2.2-2.2*Sqr(max(0,1.-sqrt(4.*g)))+(24./7.-2.2)*fg1(g, .7);
    c4fdfs=.414-.414*Sqr(max(0,1.-sqrt(2.*g)))+.086*fg3(g, 0, .5)
    +.5-.5*fg0(g, .5)+(24./7.-1.)*fg1(g, .85);
    c4rdtz=24./7. +.8*(1.-g)+(.1-24./7.)*Sqr(max(0,1.-sqrt(4.*g)));
    c4rdts=1.9-1.4*Sqr(max(0, 1-sqrt(g/.15)))-.5*fg0(g, .5)
    +(24./7.-1.9)*Sqr(g);
    
    c6fdfz=0.;
    c6fdfs=0.;
    c6rdtz=4.2-6*Sqr(max(0,1.-sqrt(4.*g)));
    c6rdts=0.;
    /*  with fw parameter */
    c3fdf=fw*c3fdfs+(1.-fw)*c3fdfz;
    c3rdt=fw*c3rdts+(1.-fw)*c3rdtz;
    c4fdf=fw*c4fdfs+(1.-fw)*c4fdfz;
    c4rdt=fw*c4rdts+(1.-fw)*c4rdtz;
    c6fdf=fw*c6fdfs+(1.-fw)*c6fdfz;
    c6rdt=fw*c6rdts+(1.-fw)*c6rdtz;
    /* calc c3 c4 c6 with fb(bfac) */
    denom=(.9*epst+.1*eps)*bb;
    if (denom==0 || modr=='n')  {capb=0; fb=0; c3=c3fdf; c4=c4fdf; c6=c6fdf;}
    else
    { 
      bparm0=(-.5*bb*((c1-2)*omm+2.*(1.+.25*c6fdf)*pdk) +.5*(c3fdf-4./3.)*bmnsmn
              +(c4fdf-2.)*bmnbmksnk) *qm*qm/denom;
      bparm1=(-.5*bb*((c1-2)*omm+2.*(1.+.25*c6rdt)*pdk) +.5*(c3rdt-4./3.)*bmnsmn
              +(c4rdt-2.)*bmnbmksnk) *qm*qm/denom;
      bplo=min(bparm0,bparm1);
      bphi=max(bparm0,bparm1);
      capb=.5*(bplo+bphi);
      for (iter=0; iter<20;iter++)
      { 
        bpold=capb; 
        fb=max(capb-.7,0.)/(max(capb-.7,0.)+1.7);
        bpnew=bparm0+fb*(bparm1-bparm0);
        bphi=max(bpnew,bpold);
        bplo=min(bpnew,bpold);
        capb=.5*(bplo+bphi);
        if (bphi<.7 || bplo>5000. || (bphi-bplo)/bphi<.001) break;
      }
      fb=max(capb-.7,0.)/(max(capb-.7,0.)+1.7); 
      c3=c3fdf*(1.-fb)+c3rdt*fb;
      c4=c4fdf*(1.-fb)+c4rdt*fb;
      c6=c6fdf*(1.-fb)+c6rdt*fb;
    }
    if (bug>0) 
    {
      beqex[mmid+iallm*lfb]=fb; 
      beqex[mmid+iallm*lc3]=c3;
      beqex[mmid+iallm*lc4]=c4; 
      beqex[mmid+iallm*lc6]=c6; 
    }
    /* include omm and c6 in c1x */
    c1xp=(1.-.5*c1)*omm;
    c1xm=0;
    if (pdk>0) c1xp -= (1+.25*c6)*pdk;
    else c1xm = -(1+.25*c6)*pdk;
    c3x=(.5*c3-2./3.);
    c4x=(.5*c4-1.);
    /* coefficients Sbij=  -beqc b[ijcen] +beqm[ij] + beqd[ij][nm] db[nmcen]
     beqc[noindfwpts], beqm[noindfwpts*6], beqd[noindfwpts*6*6] */
    tbeqc=-rhom*c1xp;    /* sign convention as in vki notes */
    Loop(n,0,6) tbeqm[n]=rhom*(c1xm*bm[n2ii[n]][n2jj[n]]
                               +c3x*sij[n2ii[n]][n2jj[n]]);
    Loop(n,0,6) Loop(k,0,6) tbeqd[n][k]=0;  
    /* material frame indifferent wall */
    frot=0;
    if (walldist>0)
    {  
      frot=walldist[wherep[In4(ip,i4dp)]]*rhom*qm/vlam[mmid];
      frot*=.0087*(.01*rt)*(.01*rt)*(.01*rt);
      if (frot>.0001) frot=1-exp(-frot);	
    }
    if (bug>0) beqex[mmid+iallm*lfr]=frot;
    Loop(n,0,6)
    { 
      i=n2ii[n]; j=n2jj[n]; nf=1+n/3;
      Loop(k,0,3)
      { 
        tbeqd[n][ij2n[i][k]] +=rhom*(c4x*sij[j][k]-wijr[j][k]-2*frot*rot[j][k]);
        tbeqd[n][ij2n[j][k]] +=rhom*(c4x*sij[i][k]-wijr[i][k]-2*frot*rot[i][k]);
        tbeqd[k][n]  += -rhom*c4x*sij[i][j]*2.*nf/3.;
        if (c5>0)   /* c5 being used as nonlinear instability parameter */
        {  tbeqd[n][ij2n[i][k]] +=rhom*.5*c5[mmid]*(wijr[j][k]+frot*rot[j][k]);
          tbeqd[n][ij2n[j][k]] +=rhom*.5*c5[mmid]*(wijr[i][k]+frot*rot[i][k]);
        }
      }
    }
      
      Loop(n,0,6) tbeqm[n] += tbeqd[n][0]*bm[0][0]+tbeqd[n][1]*bm[1][1]
      +tbeqd[n][2]*bm[2][2]+tbeqd[n][3]*bm[0][1]
      +tbeqd[n][4]*bm[0][2]+tbeqd[n][5]*bm[1][2];
      /*	printout("normal","\n j %d tbeqc %lg then tbeqm  then tbeqd \n",ip[1],tbeqc);
       Loop(i,0,6) printout("normal"," %lg",tbeqm[i]); printout("normal","\n\n");
       Loop(i,0,6) { Loop(k,0,6) printout("normal"," %lg",tbeqd[i][k]); printout("normal","\n"); }
       */
      
      /* simple c7 term not included in 6x6 matrix */
      if (c7 > 0)
      {
        Loop(n,0,6)
        {
          i=n2ii[n]; j=n2jj[n];
          Loop(k,0,3) Loop(m,0,3)
          tbeqm[n]+=.5*c7[mmid]*(bm[i][k]*bm[k][m]*sij[j][m]+bm[j][k]*bm[k][m]*sij[i][m]
                                 -2*bm[k][j]*bm[m][i]*sij[k][m]);
        }
      }
      
      Loop(k,0,8)   /* apportion by volume to equations */
      { 
        ieq=wherefw[mc[k]];
        beqc[ieq]+=tbeqc*vol[k];
        Loop(i,0,6)
        { 
          beqm[ieq+noindfwpts[0]*i]+=tbeqm[i]*vol[k];
          Loop(j,0,6) beqd[ieq+noindfwpts[0]*(i+6*j)]+=tbeqd[i][j]*vol[k];
        }
      }
    } /* ip loop */
    /* put beqc into beqm and beqd */
    Loop(i,noindfwpts[its+1],noindfwpts[ite+2])
    { 
      j=whoisfw[i];
      Loop(k,0,6) 
      { 
        beqm[i+noindfwpts[0]*k]-=beqc[i]*bij[j+iall*k];
        beqd[i+noindfwpts[0]*(k+6*k)]-=beqc[i];
      }
    }
    /* now make rhsbij = beqm + coef(0) */
    cn=(int *)need("coef_n");
    ci=(int **)need("coef_i");
    cc=(double **)need("coef_c");
    Loop(ieq,noindfwpts[its+1],noindfwpts[ite+2])
    {  
      if (cn[ieq]>0)
        Loop(i,0,cn[ieq]) Loop(j,0,6)
        rhsbij[ieq+j*noindfwpts[0]]-=bij[ci[ieq][i]+j*iall]*cc[ieq][i];
    }
    /* set cambij */	
    Loop(ieq,noindfwpts[its+1],noindfwpts[ite+2])
    { 
      cambij[ieq]=0;
      Loop(j,0,6)
      { 
        dp=0; dm=0;
        Loop(i,0,6)
        { 
          dd=beqd[ieq+noindfwpts[0]*(j+i*6)];
          if (dd<0) dm-=dd;
          else dp+=dd;
        }
        cambij[ieq]=max(cambij[ieq],max(dp,dm));
      }
    }
}