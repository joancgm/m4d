/* contains prinstress, shearijzero */
#include "global.h"

/* note these defines need available integers i,j,n,m not otherwise being used */
#define Matrixmult(a,b,c) {for (i=0;i<3;i++) for (j=0;j<3;j++) c[i][j]=a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];}
#define Copy(a,b)  {for (i=0;i<3;i++) for (j=0;j<3;j++) b[i][j]=a[i][j];} 
#define Matrixmult2(a,b,c) {for (i=0;i<3;i++) for (j=0;j<3;j++) {c[i][j]=0.; for (m=0;m<3;m++) for (n=0;n<3;n++) c[i][j] += a[i][n]*a[j][m]*b[n][m];}}
  void shearijzero(double (*re)[3], double (*tr)[3], int i, int j);
/* ---------------------------- */
double prinstress(double (*re)[3], double (*retr)[3], double (*tr)[3], double tol)                                         
/* re -   input tensor
 retr - output principle coor tensor
 tr  -  output trannsformation tensor
 tol -  input tolerance (max sum of final shears)
 */
{ 
  int i,j,n,m,iter; double err,old01,old02,old12,trn[3][3],d[3][3];
  
  for (n=0;n<3;n++) for (m=0;m<3;m++)
  {
    retr[m][n]=re[m][n];
    tr[m][n]=0; if  (m==n) tr[m][n]=1.;
  }
  
  err=abs(retr[0][1])+abs(retr[0][2])+abs(retr[1][2]);
  if (err<tol) return(err);
  for (iter=0;iter<150;iter++)
  {
    old01=abs(retr[0][1]);
    old12=abs(retr[1][2]);
    old02=abs(retr[0][2]);
    
    if (old01>= old12 && old01>= old02)    /*  0 1 transform */
    {
      shearijzero(retr,trn,0,1);
      Matrixmult(trn,tr,d); 
      Copy(d,tr);
      Matrixmult2(tr,re,retr); 
      if (abs(retr[0][1])>old01) break; /* error or roundoff */
    }
    else if (old12>=old01 && old12>=old02) /* 1 2 transform */
    { 
      shearijzero(retr,trn,1,2);
      Matrixmult(trn,tr,d); 
      Copy(d,tr);
      Matrixmult2(tr,re,retr);
      if (abs(retr[1][2])>old12) break; /* error or roundoff */
    }
    else /* 0 2 transform */
    { 
      shearijzero(retr,trn,0,2);
      Matrixmult(trn,tr,d); 
      Copy(d,tr);
      Matrixmult2(tr,re,retr);
      if (abs(retr[0][2])>old02) break; /* error or roundoff */
    }
    err=abs(retr[0][1])+abs(retr[0][2])+abs(retr[1][2]);
    /*      printout("normal",err %f tol %f \n",err,tol); */
    if (err<tol) break;
  }
  return(err);
}
/* ------------------------- */
void shearijzero(double (*re)[3], double (*tr)[3], int i, int j)
{
  int n,m; 
  double v2u2, uv2,a,b,c,cc,s,root,shear;
  
  for (n=0;n<3;n++)
  { 
    for (m=0;m<3;m++) tr[m][n]=0.; 
    tr[n][n]=1.;
  }
  if (re[i][j]==0.) return;
  
  v2u2=(re[j][j]-re[i][i])*(re[j][j]-re[i][i]);
  uv2=(re[i][j])*(re[i][j]);
  a=v2u2+4.*uv2;
  b=-v2u2-4.*uv2;
  c=uv2;
  
  root=b*b - 4.*a*c;
  cc=(-b-sqrt(root))/(2.*a);
  s=sqrt(1.-cc);
  c=sqrt(cc);
  if ((re[j][j]-re[i][i])*(cc-s*s)*re[i][j]>0) s=-s;
  shear=s*c*(re[j][j]-re[i][i])+(cc-s*s)*re[i][j];
  tr[i][i]=c;
  tr[i][j]=s;
  tr[j][i]=-s;
  tr[j][j]=c;
} 




