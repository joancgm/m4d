/* contains c_interp_gtom */
#include "global.h"

/* interpolate from grid  points to mid-points  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])

void c_interp_gtom(FILE *fpin, FILE *fprint)
{ 
  int *i4d; /* need */
  double *pg[20], *pm[20];  /* from and to */
  int n; char*name,*named; /* input parmaters */
  int i,j,iall,i4dm[4],iallm,ii[4],im;
  double ff[3]={.5,.5,.5}, *pf;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  
  n=readint(fpin);  if (n>20) n=20;
  printout("normal","interp");
  Loop(i,0,n) 
  {
    name=readname(fpin);
    pg[i]=(double *)need(name);
    named=readname(fpin);
    printout("normal",", %s to %s",name,named);
    pm[i]=(double *)find(named);
    if (pm[i]==0) 
    { 
      pm[i]=(double*)createarray(named,iallm,'d',0);
      Loop(j,0,iallm) pm[i][j]=0;
    }
  }
  printout("normal","\n");
  
  Loop(ii[3],0,i4dm[3]) Loop3(ii,0,i4dm)
  { 
    im=In4(ii,i4dm);
    Loop(j,0,n) {pf=pg[j]+ii[3]*iall; pm[j][im]= Value3d(pf,i4d,ii,ff); }
  }
}
