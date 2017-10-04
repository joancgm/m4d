/* contains c_set_xyzdouble */
#include "global.h"

/* set xyzdouble and xyzp  using geom8fx27 */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_set_xyzdouble(FILE *fpin, FILE *fprint)
{  
  int *i4d;   /* need */
  int *itrange; /* use if available */
  double *xd,*xpc; /* set xyzdouble and xyzp */
  int its,ite;  
  
  int i,j,ip[4], i4dp[4],idd[4],iallp,ialld,ilow[4],ii[4],ia[3],ipt;
  double f[27][3], x[27][3];
  int three[3]={3,3,3};
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  Loop(i,0,3) idd[i]=2*i4d[i]-1; idd[3]=i4d[3];
  ialld=Prod4(idd);
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  
  xd=(double *)createarray("xyzdouble",ialld*3,'d',0);
  Loop(i,0,ialld*3) xd[i]=0;
  xpc=(double*)createarray("xyzp",iallp*3,'d',0);
  Loop(i,0,iallp*3) xpc[i]=0;
  geom8init();
  
  Loop(ip[3],its,ite+1) 
  {
    ilow[3]=ip[3]; ii[3]=ip[3];
    Loop3(ip,1,i4d)     /* set xyzdouble */
    {
      geom8fx27(f,x,ip);
      Loop(i,0,3) ilow[i]=(ip[i]-1)*2;
      Loop3(ia,0,three)
      {
        Loop(i,0,3) ii[i]=ilow[i]+ia[i];
        ipt=In4(ii,idd);
        j=ia[0]+3*ia[1]+9*ia[2];
        Loop(i,0,3) xd[ipt+ialld*i]=x[j][i];
      }
    }
    Loop3(ip,0,i4dp)   /* set xyzp */
    {
      ipt=In4(ip,i4dp);
      Loop(i,0,3) 
      {
        ii[i]=2*ip[i]-1; 
        if (ii[i]<0) ii[i]=0; 
        else if (ii[i]>idd[i]-1) ii[i]=idd[i]-1;
      }
      j=In4(ii,idd);
      Loop(i,0,3) xpc[ipt+iallp*i]=xd[j+ialld*i];
    }
  }
}

