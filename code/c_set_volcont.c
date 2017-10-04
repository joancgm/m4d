/* contains c_set_volcont */
#include "global.h"

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_set_volcont(FILE *fpin, FILE *fprint)
{ 
  int *i4d;	 	/* needed arrays */ 
  double *volcont;    /* created  */
  
  int i4dp[4],ip[4],i,i4dm[4],mmid,im[4];
  double grad[8][3];
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  
  geomcinit();
  
  volcont=(double *)createarray("volcont",Prod4(i4dm),'d',0);
  Loop(i,0,Prod4(i4dm)) volcont[i]=0;
  
  Loop(ip[3],0,i4d[3]) Loop3(ip,1,i4d) 
  {
    Loop(i,0,3) im[i]=ip[i]-1; im[3]=ip[3];
    mmid=In4(im,i4dm);
    volcont[mmid]=geomcgradvol(ip,grad);
  }
}
