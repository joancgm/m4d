/* contains c_viscmarvheat */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* calculate companion isotropic turbulent viscosity for MARV model 
 heat transfer, vmarvheat
 see aiaapaper2010-4314  */

void c_viscmarvheat(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep;   /* need */
  double *q, *om, *rho, *prlam, *walldist, *vlam; 
  char *cltp;
  double *vmarvheat;  /* set */
  int *itrange; /* use if available */ 
  
  int i,k,i4dm[4],iallm,ip[4],i4dp[4],its,ite,ii[4],ia[3],mc[8],mmid,mwp;
  double qm,omm,rhom,cave[8],yqdnu,fac,facpr;
  int two[3]={2,2,2};
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  
  rho=(double *)need("rho");
  q=(double *)need("qturb");
  om=(double *)need("omturb");
  prlam=(double *)need("prlam");
  vlam=(double *)need("vlam");
  walldist=(double *)need("walldist");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1; }
  
  geomcinit();
  vmarvheat=(double *)createarray("vmarvheat",iallm,'d',0);
  
  /* loop over contiuity control volumes */
  Loop(ip[3],its,ite+1) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));  /* viscosity index */      vmarvheat[mmid]=0;   /* initialize as zero */
    mwp=wherep[In4(ip,i4dp)];
    if (mwp<0) continue;  /* not valid c.v. leave as zero */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }  
    geomcvave(cave,ip);   /* real cont. c.v. ave */
    qm=0; omm=0; rhom=0;
    Loop(k,0,8)
    { 
      rhom+=cave[k]*rho[mc[k]];
      /* omm+=cave[k]/om[mc[k]]; */ /* assume 1/om varies linearly */
      omm+=cave[k]*log(om[mc[k]]);
      qm+=cave[k]*q[mc[k]];
    }
    /* omm=1./omm; */
    omm=exp(omm); 
    yqdnu=walldist[mwp]*qm/(vlam[mmid]/rhom);
    if (yqdnu<=0) continue;  /*  leave as zero */
    fac=.00018*yqdnu*yqdnu;
    facpr=fac*pow(prlam[0],0.7);
    if (fac>1.e-5) fac=1-exp(-fac);
    if (facpr>1.e-5) facpr=1-exp(-facpr);
    vmarvheat[mmid]=.1*sqrt(fac*facpr)*rhom*qm*qm/omm;
  }
}
