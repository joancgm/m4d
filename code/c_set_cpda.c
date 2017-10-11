/* contains c_set_contar, c_set_cpda */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

/* ------------------------- */
/* determine aspect ratios of continuity control volumes */
void c_set_contar(FILE *fpin, FILE *fprint)
{ 
  int *noindppts,*whoisp,*i4d; double *x[3]; /* needed arrays */
  double *contar; /* create */
  int i,j,k,n,idpc[4],ip[4],iall,iallp,iadd[3],ia[3],ipg,igg;
  double dx[3][3],da[3],aa[3],rmin,rmax;
  int nopts;
  
  /*printout("normal","enter contar\n"); */
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  iadd[0]=1; iadd[1]=i4d[0]; iadd[2]=i4d[0]*i4d[1];
  
  Loop(n,0,3) idpc[n]=i4d[n]+1; idpc[3]=i4d[3];
  iallp=Prod4(idpc);
  noindppts=(int *)need("noindppts");
  nopts=noindppts[0];
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[1]+iall;
  whoisp=(int *)need("whoisp");
  contar=(double *)createarray("contar",nopts*3,'d',0);
  
  Loop(i,0,nopts)
  {
    iexpand(whoisp[i],idpc,ip);
    Loop(j,0,3) contar[i+nopts*j]=0;
    if (ip[0]==0 || ip[0]==i4d[0] || ip[1]==0 || ip[1]==i4d[1] 
        || ip[2]==0 || ip[2]==i4d[2]) continue;
    ipg=In4(ip,i4d);
    /*printout("normal","i %d ip %d %d %d %d ipg %d iall %d\n",i,ip[0],ip[1],ip[2],ip[2],ipg,iall); */
    Loop(j,0,3) Loop(k,0,3) dx[j][k]=0;
    Loop(ia[0],0,2) Loop(ia[1],0,2) Loop(ia[2],0,2)
    { 
      igg=ipg-ia[0]*iadd[0]-ia[1]*iadd[1]-ia[2]*iadd[2];
      /*printout("normal"," igg %d  x %lg %lg %lg\n",igg,x[0][igg],x[1][igg],x[2][igg]); */
      Loop(j,0,3) Loop(k,0,3)
      dx[j][k]+=x[k][igg]*(2*ia[j]-1);
    }
    /*Loop(j,0,3) Loop(k,0,3) printout("normal"," %lg ",dx[j][k]); printout("normal","\n"); */
    Loop(j,0,3) aa[j]=Sqr(dx[j]);
    rmin=min(aa[0],min(aa[1],aa[2]));
    if (rmin==0) /* use length */
    { 
      Loop(j,0,3) if (aa[j]!=0) aa[j]=sqrt(1/aa[j]);  /* 1/dx */
    }
    else /* use normal width */
    { 
      Loop(j,0,3) 
      {
        Cross(dx[j],dx[(j+1)%3],da);
        /*printout("normal","j %d da %lg %lg %lg\n",(j+2)%3,da[0],da[1],da[2]);*/
        aa[(j+2)%3]=sqrt(Sqr(da));
      }
    }
    rmax=max(aa[0],max(aa[1],aa[2])); /* max of 1/dx(non-zero) */
    if (rmax>0) Loop(j,0,3) contar[i+nopts*j]=aa[j]/rmax;
    /*printout("normal","aa %lg %lg %lg rmax %lg contar %lg %lg %lg\n",
     aa[0],aa[1],aa[2],rmax,contar[i],contar[i+nopts],contar[i+2*nopts]); */
  }
  printout("normal","set contar\n");
}
/* ------------------------- */
/*  coefficients pressure term int pdAi  for momentum equations in terms of P at pc points,
 set up coefficients, cpda   for specified time indices */
void c_set_cpda(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep,*noindppts;    /* needed arrays */
  double  *xpc; char *clt,*csym; 
  int *itrange; double *roundoff; /* use if available */
  double *cpda;  /* create or modify */
  double *contar; /* find or set */
  double cenfac[3],ar1,ar2; /* input parameters */
  int nomit; char omit[20]; 
  
  int it,iall,idd[4],ialld,idpc[4],iallp,nyclear;
  int i,j,k,n,ip[4],L,L2,L3,ig[4],id[4],ja,jb,jc[3],jg,ipa[4],ia[3],ippt;
  int its,ite;
  double xyzc[3],area[3],xyzpc[2][2][3],fpc[2],ff,ffa,ffb;
  double ar,cenf[3]={0,0,0};
  double tol=1.e-11;
  int icbug=0;  
  double areak[20];
  Loop(i,0,20) areak[i]=0;
  
  cenfac[0]=readdouble(fpin);   /* read input */
  ar1=readdouble(fpin);
  cenfac[1]=readdouble(fpin);
  ar2=readdouble(fpin);
  cenfac[2]=readdouble(fpin);
  Loop(i,0,3) cenfac[i]=max(0,min(1,cenfac[i]));
  if (ar2<ar1)
  {
    ar=ar1; ar1=ar2; ar2=ar; 
    ar=cenfac[1]; cenfac[2]=cenfac[1]; cenfac[1]=ar;
  }
  printout("normal","ar=0 cenfac=%lg, linear: ar=%lg cenfac=%lg to ar=%lg cenfac=%lg\n",
          cenfac[0],ar1,cenfac[1],ar2,cenfac[2]);
  nomit=readint(fpin);
  printout("normal"," omit %d pt types:",nomit);
  Loop(i,0,nomit) 
  {
    omit[i]=read1charname(fpin); 
    printout("normal"," %c",omit[i]); 
  }
  printout("normal","\n");
  
  contar=(double *)find("contar");   /* make sure contar exists */
  if (contar==0) 
  {
    c_set_contar(fpin,fprint);
    contar=(double *)need("contar");
    icbug+=1;
  }
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  noindppts=(int *)need("noindppts");
  xpc=(double *)need("xyzp");
  clt=(char *)need("clt");
  csym=(char *)need("csym");
  cpda=(double *)find("cpda"); 
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  geom8init();  
  
  Loop(n,0,3) { idd[n]=2*i4d[n]-1;  idpc[n]=i4d[n]+1;}
  idd[3]=i4d[3]; idpc[3]=i4d[3];
  iall=Prod4(i4d);
  ialld=Prod4(idd);
  iallp=Prod4(idpc); 
  /* printout("normal"," iall: %d %d %d\n", iall,ialld,iallp); */
  if (cpda==0)  
  { 
    cpda=(double *)createarray("cpda",8*3*iall,'d',0);
    Loop(n,0,8*3*iall) cpda[n]=0;
  }
  else  /* clear ones to  be reset */
  {
    Loop(i,0,3) ip[i]=0; 
    ip[3]=its; 	 ja=24*In4(ip,i4d);
    ip[3]=ite+1;  jb=24*In4(ip,i4d);
    Loop(i,ja,jb) cpda[i]=0;
  }
  
  for (it=its;it<=ite;it++)   /* time index */
  { 
    ip[3]=it;  ig[3]=it;  id[3]=it; ipa[3]=it;
    Loop(L,0,3)   /* each direction */
    {
      L2=(L+1)%3; L3=(L+2)%3; ia[L]=0;
      Loop(ip[L],0,i4d[L]+1)
      Loop(ip[L2],1,i4d[L2]) Loop(ip[L3],1,i4d[L3]) /* each pc volume */
      {
        ippt=In4(ip,idpc);
        /* active pc points or end boundaries  only */
        if (wherep[ippt] >=0 || ip[L]==0 || ip[L]==i4d[L]) 
        {
          Loop(i,0,3)
          {
            if (i==L) continue;
            ar=0; 
            if (wherep[ippt]>=0) ar=contar[wherep[ippt]+i*noindppts[0]];
            if (ar==0) cenf[i]=cenfac[0];
            else if (ar<=ar1) cenf[i]=cenfac[1];
            else if (ar>=ar2) cenf[i]=cenfac[2];
            else cenf[i]=cenfac[1]+(ar-ar1)*(cenfac[2]-cenfac[1])/(ar2-ar1);
            if (icbug>1 &&ar>0 && ar<ar2) 
              printout("normal","ip %d %d %d L %d i %d ar %lg cenf %lg\n",
                     ip[0],ip[1],ip[2],L,i,ar,cenf[i]);
          }
          Loop(ia[L2],0,2) Loop(ia[L3],0,2) /* each quarter area */
          {
            id[L]=max(0,min(idd[L]-1,2*ip[L]-1));
            id[L2]=2*ip[L2]-2+ia[L2];
            id[L3]=2*ip[L3]-2+ia[L3];
            geom8areamid(area,xyzc,id,L);  /* get area and midpoint */
            ipa[L]=ip[L];
            Loop(ja,0,2) Loop(jb,0,2)  /* get xyzpc */
            { 
              ipa[L2]=ip[L2]+ja+ia[L2]-1;
              ipa[L3]=ip[L3]+jb+ia[L3]-1;
              k=In4(ipa,idpc);
              Loop(n,0,3) xyzpc[ja][jb][n]=xpc[k+iallp*n];
            }
            fpc[0]=.5; fpc[1]=.5;
            k=  find2dinterp(xyzc,xyzpc,fpc,1,2);  /* returns fpc */
            Loop(i,0,3) ig[i]=ip[i]+ia[i]-1;
            Loop(ja,0,2) Loop(jb,0,2)  /* add contributions to cpda */
            { 
              jc[L2]=ja; jc[L3]=jb;
              /* direction biased ar based pressure for area */
              ffa=(ja*fpc[0]+(1-ja)*(1-fpc[0]));
              if (ja+ia[L2]-1==0) ffa=cenf[L2]+(1-cenf[L2])*ffa;
              else ffa=(1-cenf[L2])*ffa;
              ffb=(jb*fpc[1]+(1-jb)*(1-fpc[1]));
              if (jb+ia[L3]-1==0) ffb=cenf[L3]+(1-cenf[L3])*ffb;
              else ffb=(1-cenf[L3])*ffb;
              ff=ffa*ffb;
              /* end try */
              ig[L]=ip[L]-1;
              if (ig[L]>=0)
              { 
                jc[L]=1;
                j=jc[0]+2*(jc[1]+2*jc[2]);
                jg=In4(ig,i4d);
                Loop(n,0,3) 
                cpda[j+8*n+24*jg]+=area[n]*ff;
              }
              ig[L]=ip[L];
              if (ig[L]<i4d[L])
              {
                jc[L]=0;
                j=jc[0]+2*(jc[1]+2*jc[2]);
                jg=In4(ig,i4d);
                Loop(n,0,3) 
                cpda[j+8*n+24*jg]-=area[n]*ff;
              }
            }  /* ja jb */
          }   /* ia loop */
        } /* match */
      } /* each pc volume */
    } /* Loop L */
  } /* loop it */
  jb=0; n=0;
  /* clear values for wall, solid, other specified points and as appropriate on sym planes */
  Loop(ig[3],its,ite+1) Loop3(ig,0,i4d)
  {
    i=In4(ig,i4d);
    nyclear=0;
    if (clt[i]=='w' || clt[i]=='s') nyclear=1;
    else Loop(j,0,nomit) if (clt[i]==omit[j]) nyclear=1;
    if (nyclear==1)
    { 
      jb++;
      Loop(ja,0,24) cpda[ja+24*i]=0;
      continue;
    }
    /* check for sym planes */
    Loop(j,0,6) if (csym[j]!='n')
    {
      if (j==0 && ig[0]!=0) continue;
      if (j==1 && ig[0]!=i4d[0]-1) continue;
      if (j==2 && ig[1]!=0) continue;
      if (j==3 && ig[1]!=i4d[1]-1) continue;
      if (j==4 && ig[2]!=0) continue;
      if (j==5 && ig[2]!=i4d[2]-1) continue; 
      L=-1;
      if (csym[j]=='x') L=0;
      else if (csym[j]=='y') L=1;
      else if (csym[j]=='z') L=2;
      
      if (L>=0) 
      {
        n++;
        Loop(k,0,8) cpda[k+8*L+24*i]=0;
      }
    }
    
  }
  
  if (jb>0) printout("normal"," cpda cleared at %d points based on type\n",jb);
  if (n>0) printout("normal","   components cleared at %d points based on csym\n",n);
  /* clear roundoff */
  for (i=0;i<8*3*iall;i+=8)
  { 
    ff=0;
    for (j=0;j<8;j++) ff=max(ff,abs(cpda[i+j]));
    ff*=tol;
    for (j=0;j<8;j++)
      if (abs(cpda[i+j])<ff) cpda[i+j]=0;
  }
}
