/* contains c_contdu */
#include "global.h"

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])

/* calc  dU from dp and update U */
void c_contdu(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*match,*wherefw,*cpflop_n, **cpflop_i; /* needed arrays */
  double *cam,*cplus,**cpflop_c,*dp;
  int *itrange; /* use if available */
  double *u[3], *du[3];     /* update */
  int ii[4],ipt,iw,i,j,n;
  double psum;
  int its,ite;
  
  i4d=(int *)need("idim4d");
  match=(int *)need("match");
  wherefw=(int *)need("wherefw"); 
  cpflop_n=(int*)need("cpflop_n"); 
  cpflop_i=(int**)need("cpflop_i"); 
  cpflop_c=(double**)need("cpflop_c"); 
  cam=(double *)need("cam");
  cplus=(double *)need("cplus");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  du[0]=(double *)find("dU1");
  du[1]=(double *)find("dU2");
  du[2]=(double *)find("dU3");
  dp=(double *)need("dp");
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  
  Loop(ii[3],its,ite+1) Loop3(ii,0,i4d)  /* do active points first */
  { 
    ipt=In4(ii,i4d); iw=wherefw[ipt];
    if (iw>=0)
    { 
      n=cpflop_n[iw];
      if (n>0 && cplus[iw]>0) 
        Loop(j,0,3)
      { 
        psum=0;
        Loop(i,0,cpflop_n[iw]) psum+=cpflop_c[iw][i+n*j]*dp[cpflop_i[iw][i]];
        u[j][ipt]-=psum/(cam[iw]+cplus[iw]);
        if (du[j]>0) du[j][ipt]-=psum/(cam[iw]+cplus[iw]);
      }
    }
  }
  Loop(ii[3],its,ite) Loop3(ii,0,i4d)  /* fill in passive points */
  { 
    ipt=In4(ii,i4d);
    Loop(i,0,3) u[i][ipt]=u[i][match[ipt]];
  }
}
