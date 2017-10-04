/* contains c_coefdt */
#include "global.h"
/* calculate coefficient for time term  */

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_coefdt(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep,*noindfwpts,*wherefw,*match,*whoisfw,*nocoefs;     /* needed arrays */
  double *dt,*rho;
  char *cltp;
  int *itrange; double *roundoff; /* use if available */
  int *coef_n, **coef_i; double **coef_c;  /* modify*/
  char *namedt; int jdt; char nymidbias;  /* input parameters */
  
  int i,j,its,ite,ip[4],ig[4],jv,ia[3],iw[8],iptc[8],i4dp[4],i4dm[4],iallm,nylocal;
  double vol[8],fmid[8][3];
  double tol=1.e-8;
  int two[3]={2,2,2};
  double ffc[8],rhomvdt,ff[3],dtc;
  int bug=0;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  
  namedt=readname(fpin);
  jdt=readint(fpin);
  nymidbias=read1charname(fpin);
  
  dt=(double *)need(namedt);
  nylocal=0;
  if (arraysize(namedt)==iallm) nylocal=1;
  printout("normal"," dt: %s",namedt);
  if (nylocal==0) printout("normal"," = %lg,",dt[0]);
  else printout("normal",", local values,");
  printout("normal"," at coef %d, nymidbias %c",jdt,nymidbias);
  
  nocoefs=(int *)need("nocoefs");
  if (jdt>nocoefs[0]-1)
  { 
    printout("error c_coefdt"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
    exitm4d(0);
  }
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  rho=(double *)need("rho"); 
  coef_n=(int *)need("coef_n");
  coef_i=(int **)need("coef_i");
  coef_c=(double **)need("coef_c");
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  noindfwpts=(int *)need("noindfwpts");
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  wherefw=(int *)need("wherefw");
  whoisfw=(int *)need("whoisfw");
  match=(int *)need("match");
  geom8init();  
  
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])  /* clear parts of array to be redone */
  if (coef_n[i]>0) {Loop(j,0,coef_n[i]) coef_c[i][j+coef_n[i]*jdt]=0; }
  
  Loop(ip[3],its,ite+1) /* loop over time steps */
  { 
    ig[3]=ip[3];
    Loop3(ip,1,i4d) /* loop over all possible cont c.v. */
    { 
      if (wherep[In4(ip,i4dp)]<0) continue; /* skip if not valid cont c.v. */
      if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
      dtc=dt[0];
      if (nylocal==1) dtc=dt[ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]))];
      if (dtc==0) continue;
      geom8volfmid(vol,fmid,ip);
      Loop3(ia,0,two) /* each corner */
      { 
        jv=ia[0]+2*ia[1]+4*ia[2];
        Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
        i=In4(ig,i4d);
        iptc[jv]=match[i];
        iw[jv]=wherefw[iptc[jv]];
        if (bug>0) printout("normal"," jv i iptc iw %d %d %d %d \n",jv,i,iptc[jv],iw[jv]);
      }
      Loop3(ia,0,two)  /* each corner volume */
      {
        jv=ia[0]+2*ia[1]+4*ia[2];
        Loop(i,0,3) ff[i]=fmid[jv][i];
        if (bug>0)printout("normal","jv %d  iw %d  ff %lg %lg %lg\n",jv,iw[jv],ff[0],ff[1],ff[2]);
        if (iw[jv]<0) continue;
        ffc[0]=(1.-ff[0])*(1.-ff[1])*(1.-ff[2]);
        ffc[1]=ff[0]*(1.-ff[1])*(1.-ff[2]);
        ffc[2]=(1.-ff[0])*ff[1]*(1.-ff[2]);
        ffc[3]=ff[0]*ff[1]*(1.-ff[2]);
        ffc[4]=(1.-ff[0])*(1.-ff[1])*ff[2];
        ffc[5]=ff[0]*(1.-ff[1])*ff[2];
        ffc[6]=(1.-ff[0])*ff[1]*ff[2];
        ffc[7]=ff[0]*ff[1]*ff[2];
        rhomvdt=0;
        Loop(i,0,8) rhomvdt+=ffc[i]*rho[iptc[i]];
        rhomvdt*=vol[jv]/dtc;
        Loop(i,0,8) ffc[i]*=rhomvdt;
        coefcadd(iw[jv],coef_n,coef_i,coef_c,nocoefs[0],jdt,8,iptc,ffc); 
      }
    } /* ip */
  } /* ip[3] */
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])    /* clean up redundent coefficients */
  {
    coefcombine(i,coef_n,coef_i,coef_c,nocoefs[0],tol);
    coeforder(whoisfw[i],i,coef_n,coef_i,coef_c,nocoefs[0]);  
  }  /* add opposite fold to reg coefdt, simple without repeat or sleeping fold */
               if (nymidbias=='y')
               { 
                 int icen[4],ilow[4],ihi[4],jlow,jhi,klow,khi,k;
                 double dc;
                 Loop(i,noindfwpts[its+1],noindfwpts[ite+2])
                 {  
                   if (coef_n[i]<=2) continue;
                   iexpand(whoisfw[i],i4d,icen);
                   Loop3(ia,0,two)
                   { 
                     if (ia[0]==0 && ia[1]==0 && ia[2]==0) continue;
                     Loop(k,0,4) {ilow[k]=icen[k]; ihi[k]=icen[k]; }
                     Loop(k,0,3) {ihi[k]+=ia[k]; ilow[k]-=ia[k]; }
                     jlow=In4(ilow,i4d);
                     klow=-1;
                     Loop(k,0,coef_n[i])	if (jlow==coef_i[i][k]) {klow=k; break; }
                     if (klow<0) continue;
                     jhi=In4(ihi,i4d);
                     khi=-1;
                     Loop(k,0,coef_n[i])	if (jhi==coef_i[i][k]) {khi=k; break; }
                     if (khi<0) continue;
                     dc=min(coef_c[i][klow+coef_n[i]*jdt],coef_c[i][khi+coef_n[i]*jdt]);
                     if (dc<=0) continue;
                     coef_c[i][klow+coef_n[i]*jdt]-=dc;
                     coef_c[i][khi+coef_n[i]*jdt]-=dc;
                     coef_c[i][coef_n[i]*jdt]+=dc+dc;
                   }
                 }
               }
}

