/* contains c_ppreset, ppreset */
#include "global.h"
#include "psleepcalc.h"
/* reset p for moved grid */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
void ppreset(double *pp, double *xpc, double *xpcold, FILE *fprint);

void c_ppreset(FILE *fpin, FILE *fprint)
{
  double *xpc,*xpcold,*pp; 
  xpc=(double *)need("xyzp");
  xpcold=(double *)need("xyzpold");
  pp=(double *)need("pp");
  ppreset(pp,xpc,xpcold,fprint);
}
/* --------------------------- */
void ppreset(double *pp, double *xpc, double *xpcold, FILE *fprint)
{ 
  int *i4d,*matchpc;   /* need */
  double *cpsleep;
  int *itrange; /* use if available */
  int its,ite; 
  int i,ip[4],i4dp[4],iallp,ilo[3],ihi[3],iapp,ipt;
  double *ptemp,f[3],*x[3],xx[3];
  double tol=.001;
  int iter;
  
  i4d=(int *)need("idim4d");
  matchpc=(int *)need("matchpc");
  cpsleep=(double *)need("cpsleep");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1]; }
  else {its=0; ite=i4d[3]-1;}
  
  ptemp=(double *)tmalloca(iallp,'d');   /* copy pp before update */
  Loop(i,0,iallp) ptemp[i]=pp[i];
  
  Loop(ip[3],its,ite+1)   /* interpolation update */
  {
    iapp=ip[3]*i4dp[0]*i4dp[1]*i4dp[2]; 
    x[0]=xpcold+iapp; 
    x[1]=x[0]+iallp; 
    x[2]=x[1]+iallp;
    Loop3(ip,0,i4dp) 
    { 
      ipt=In4(ip,i4dp);
      /*  printout("normal","ipt %d ptemp %lg pp %lg\n",ipt,ptemp[ipt],pp[ipt]); */
      if (matchpc[ipt]>=0) continue;  /* these points done by psleepcalc */
      Loop(i,0,3)
      { 
        if (ip[i]==0) {ilo[i]=0; ihi[i]=0; }
        else if (ip[i]==i4dp[i]-1) {ilo[i]=i4dp[i]-1; ihi[i]=i4dp[i]-1; }
        else 
        {
          ilo[i]=max(ip[i]-1,1);
          ihi[i]=min(ip[i]+1,i4dp[i]-2);
        }
        f[i]=ip[i]; 
        xx[i]=xpc[ipt+i*iallp];
      }
      /*     printout("normal","ip %d %d %d %d ilo %d %d %d  ihi %d %d %d  f %lg %lg %lg\n",
       ip[0],ip[1],ip[2],ip[3],ilo[0],ilo[1],ilo[2],ihi[0],ihi[1],ihi[2], f[0],f[1],f[2]);  */
      iter=findex3d(xx,x,i4dp,ilo,ihi,tol,f);
      interp3d(ptemp+iapp,i4dp,f,pp+ipt);
      /*	   printout("normal"," iter %d f %lg %lg %lg pold %lg pnew %lg\n",
       iter,f[0],f[1],f[2],ptemp[ipt],pp[ipt]);  */
    }
  }
  tmalloca(-1,'i');
  psleepcalc(fprint,pp,cpsleep);
}


