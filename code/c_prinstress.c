/* contains c_prinstress */
#include "global.h"

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Loop(n,a,b) for (n=a;n<b;n++)

void c_prinstress(FILE *fpin, FILE *fprint)
/*  calc principal arrays, Reynolds stresses, strain rate and vorticity
 which can be used by a.x11pict */
{ 
  double *uddxi[3], *bij, *q;    /* input create U1(2,3)ddxi with gradprop */ 
  int *i4d;
  
  double *repr, *stpr, *vvec;    /* output */
  
  int i,iall,iallm,i4dm[4],k,n,io[3],ii[4];
  double re[3][3],retr[3][3],tr[3][3],err,tol,du[3][3];
  double prinstress(double (*re)[3], double (*retr)[3], double (*tr)[3], double tol);
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallm=Prod4(i4dm);
  uddxi[0]=(double*)find("U1ddxi");
  uddxi[1]=(double*)find("U2ddxi");
  uddxi[2]=(double*)find("U3ddxi");
  q=(double*)find("qturb");
  bij=(double*)find("bij");
  
  /* Reynolds stress tensor principal vectors */
  if (q>0 && bij>0)
  {
    printout("normal"," setting Repr");
    repr=(double*)createarray("Repr",9*iall,'d',0);
    Loop(i,0,9*iall) repr[i]=0;
    ii[3]=0; /* do for first can later read if need to */
    Loop3(ii,0,i4d)
    { 
      i=In4(ii,i4d);
      if (q[i]==0) continue;
      
      for (n=0;n<3;n++) re[n][n]=2.*q[i]*q[i]*(bij[i+iall*n]+1./3.);
      re[0][1]=2.*q[i]*q[i]*bij[i+iall*3];       re[1][0] = re[0][1];
      re[0][2]=2.*q[i]*q[i]*bij[i+iall*4];       re[2][0] = re[0][2];
      re[1][2]=2.*q[i]*q[i]*bij[i+iall*5];       re[2][1] = re[1][2];
      err=prinstress(re,retr,tr,1.e-6*q[i]*q[i]);
      io[0]=0; 
      if (retr[1][1]> retr[0][0]) io[0]=1; 
      if (retr[2][2]>retr[io[0]][io[0]]) io[0]=2;
      io[2]=2;
      if (retr[1][1]< retr[2][2]) io[2]=1; 
      if (retr[0][0]<retr[io[2]][io[2]]) io[2]=0;
      io[1]=3-io[0]-io[2];
      /* tr = tr[prinvector][xyzcomponent] */
      for (n=0;n<3;n++) for (k=0;k<3;k++) repr[k+3*n+9*i]=retr[io[n]][io[n]]*tr[io[n]][k];
      if (err>1.e-5*q[i]*q[i]) 
      { 
        printout("normal","i,j,k,q %d %d %d %f  io %d %d %d\n",ii[0],ii[1],ii[2],q[i],io[0],io[1],io[2]);
        printout("normal","re");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",re[n][k]); printout("normal","\n");
        printout("normal","retr");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",retr[n][k]); printout("normal","\n");
        printout("normal","tr");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",tr[n][k]); printout("normal","\n");
      }
    }
  }
  /* Strain rate tensor principal vectors, and vorticity*/
  if (uddxi[0]>0 && uddxi[1]>0 && uddxi[2]>0)
  { 
    printout("normal","  setting Stpr, VVec");
    stpr=(double*)createarray("Stpr",9*iallm,'d',0);
    vvec=(double*)createarray("Vvec",3*iallm,'d',0);
    Loop(i,0,9*iallm)  stpr[i]=0;
    Loop(i,0,3*iallm)  vvec[i]=0;
    ii[3]=0;
    Loop3(ii,0,i4dm)
    { 
      tol=0; 
      i=In4(ii,i4dm);
      Loop(n,0,3) Loop(k,0,3) 
      {
        du[k][n]=uddxi[k][i+n*iallm]; /* dUk dxn*/
        tol=max(tol,abs(du[k][n]));
      }
      if (tol==0) continue;
      re[0][0]=du[0][0]; re[1][1]=du[1][1]; re[2][2]=du[2][2];
      re[0][1]=.5*(du[0][1]+du[1][0]);  re[1][0]=re[0][1];
      re[0][2]=.5*(du[0][2]+du[2][0]);  re[2][0]=re[0][2];
      re[1][2]=.5*(du[1][2]+du[2][1]);  re[2][1]=re[1][2];
      err=prinstress(re,retr,tr,1.e-6*tol);
      io[0]=0; 
      if (retr[1][1]> retr[0][0]) io[0]=1; 
      if (retr[2][2]>retr[io[0]][io[0]]) io[0]=2;
      io[2]=2;
      if (retr[1][1]< retr[2][2]) io[2]=1; 
      if (retr[0][0]<retr[io[2]][io[2]]) io[2]=0;
      io[1]=3-io[0]-io[2];
      /* tr = tr[prinvector][xyzcomponent] */
      for (n=0;n<3;n++) for (k=0;k<3;k++) stpr[k+3*n+9*i]=retr[io[n]][io[n]]*tr[io[n]][k];
      if (err>1.e-5*tol) 
      {
        printout("normal","i,j,k,tol %d %d %d %f  io %d %d %d\n",ii[0],ii[1],ii[2],tol,io[0],io[1],io[2]);
        printout("normal","re");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",re[n][k]); printout("normal","\n");
        printout("normal","retr");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",retr[n][k]); printout("normal","\n");
        printout("normal","tr");  for (n=0;n<3;n++) for (k=0;k<3;k++) printout("normal"," %g",tr[n][k]); printout("normal","\n");
      }
      /* vorticity */
      vvec[3*i]=du[2][1]-du[1][2];
      vvec[1+3*i]=du[0][2]-du[2][0];
      vvec[2+3*i]=du[1][0]-du[0][1];
    }
  }
  printout("normal","\n");
}

