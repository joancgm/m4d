/* contains c_copy_ptom */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_copy_ptom(FILE *fpin, FILE *fprint)   
/*  copy p-point to mid-point for arrays of type double */
{ 
  char *namefrom, *nameto; 
  int *i4d;
  int i4dp[4],i4dm[4],i,ip[4],im[4],iallm,iallp;
  double *vf,*vto;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iallp=Prod4(i4dp);
  
  namefrom=readname(fpin);
  nameto=readname(fpin);
  printout("normal","copy %s to %s\n",namefrom,nameto);
  i=arraysize(namefrom);
  if (i != iallp) 
  {
    printout("error c_copy_ptom"," array %s has size %d not %d, cannot copy\n",namefrom,i,iallp); 
    exitm4d(0); 
  }
  vf=(double *)need(namefrom);
  vto=(double *)createarray(nameto,iallm,'d',0);
  
  Loop(im[3],0,i4dm[3]) Loop3(im,0,i4dm)
  {  
    Loop(i,0,3) ip[i]=im[i]+1; ip[3]=im[3];
    vto[In4(im,i4dm)]=vf[In4(ip,i4dp)];
  }
}
