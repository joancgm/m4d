/* contains c_interp_gtop */
#include "global.h"
#include "psleepcalc.h"

/* interpolate from grid  points to p  points  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])

void c_interp_gtop(FILE *fpin, FILE *fprint)
{ 
  int *i4d; /* need */
  double *cvdc, *cpsleep; /* need use for interp */
  double *pg[20], *pp[20];  /* from and to */
  int n; char*name,*named; /* input parmaters */
  int i,j,iall,i4dp[4],iallp,ii[4],iic[4],ic,ip;
  double ff[3], *pf;
  
  i4d=(int *)need("idim4d");
  cvdc=(double *)need("cvdc");
  cpsleep=(double *)need("cpsleep");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  
  n=readint(fpin);  if (n>20) n=20;
  printout("normal","interp");
  Loop(i,0,n) 
  {
    name=readname(fpin);
    pg[i]=(double *)need(name);
    named=readname(fpin);
    printout("normal",", %s to %s",name,named);
    pp[i]=(double *)find(named);
    if (pp[i]==0) 
    { 
      pp[i]=(double*)createarray(named,iallp,'d',0);
      Loop(j,0,iallp) pp[i][j]=0;
    }
  }
  printout("normal","\n");
  
  Loop(ii[3],0,i4dp[3]) Loop3(ii,0,i4dp)
  { 
    iic[3]=0;
    Loop(i,0,3) ff[i]=-1;
    Loop(i,0,3) 
    {
      if (ii[i]==0) {iic[i]=0; ff[i]=0;}
		else if (ii[i]==i4d[i]) {iic[i]=i4d[i]-2; ff[i]=1;}
		else iic[i]=ii[i]-1;
    }
    ic=In4(iic,i4d);
    ip=In4(ii,i4dp);
    Loop(i,0,3) if (ff[i]<0) ff[i]=cvdc[ic+iall*i];
    Loop(j,0,n) {pf=pg[j]+ii[3]*iall; pp[j][ip]= Value3d(pf,i4d,iic,ff); }
  }
  Loop(j,0,n) psleepcalc(fprint,pp[j],cpsleep);
}
