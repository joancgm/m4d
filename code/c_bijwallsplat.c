/* contains c_bijwallsplat */
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Loop(n,a,b) for (n=a;n<b;n++)

void c_bijwallsplat(FILE *fpin, FILE *fprint)  

/* set bij at walls consistent with near wall and wall reflection */
{
  int   *i4d,*noindwpts,*whoisw,*wnear,**whoelse;     /* need */
  double  *wnorm; 
  int *itrange; /* use if available */
  double *bij; /* modify */ 
  
  int its,ite;
  int iall,i,j,k,m,n,nn,inw;
  double up[3],*wn,umag,dbn,ub[3],bnpb[3][3],bijw[3][3],*wt[3];
  int bug=0;
  
  i4d=(int *)need("idim4d");             /* find input values and arrays */
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  itrange=(int *)find("itrange");
  if (itrange>0) {its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1;}
  wnear=(int *)need("wnear");
  noindwpts=(int *)need("noindwpts");
  whoisw=(int *)need("whoisw");
  whoelse=(int **)need("whoelse");
  wnorm=(double *)need("wnorm");
  bij=(double *)need("bij");
  wt[1]=up; wt[2]=ub;
  Loop(i,noindwpts[its+1],noindwpts[ite+2])
  { 
    m=whoisw[i];
    inw=wnear[i];
    if (bug>0) printout("normal","m,inw %d %d\n",m,inw); 
    wn=wnorm+3*i; 
    wt[0]=wn;
    up[0]=1; up[1]=0; up[2]=0;
    if (bug>0) printout("normal"," un %lg %lg %lg up %lg %lg %lg  Dot %lg\n",
                      wn[0],wn[1],wn[2],up[0],up[1],up[2],Dot(wn,up)); 
    if (abs(Dot(wn,up))>.999) {up[0]=0; up[1]=1;}
    Cross(wn,up,ub); 
    umag=sqrt(Dot(ub,ub));
    if (bug>0) printout("normal"," up %lg %lg %lg ub %lg %lg %lg umag %lg\n",
                      up[0],up[1],up[2],ub[0],ub[1],ub[2],umag);
    Loop(n,0,3) ub[n]/=umag;
    Cross(ub,wn,up);
    if (bug>0) printout("normal"," un %lg %lg %lg up %lg %lg %lg ub %lg %lg %lg\n",
                      wn[0],wn[1],wn[2],up[0],up[1],up[2],ub[0],ub[1],ub[2]); 
    Loop(n,0,3) bijw[n][n]=bij[inw+iall*n];
    bijw[0][1]=bij[inw+iall*3]; bijw[1][0]=bijw[0][1];
    bijw[0][2]=bij[inw+iall*4]; bijw[2][0]=bijw[0][2];
    bijw[2][1]=bij[inw+iall*5]; bijw[1][2]=bijw[2][1];
    if (bug>0)  printout("normal","bijw %lg %lg %lg %lg %lg %lg\n",
                       bijw[0][0],bijw[1][1],bijw[2][2],bijw[0][1],bijw[0][2],bijw[1][2]); 
    Loop(j,0,3) Loop(k,0,3)
    { 
      bnpb[j][k]=0;
      Loop(n,0,3) Loop(nn,0,3) bnpb[j][k]+=wt[j][n]*wt[k][nn]*bijw[n][nn];
    }
    if (bug>0)  printout("normal","bnpb %lg %lg %lg %lg %lg %lg\n",
                       bnpb[0][0],bnpb[1][1],bnpb[2][2],bnpb[0][1],bnpb[0][2],bnpb[1][2]); 
    dbn=-.33333-bnpb[0][0];
    if (dbn<0)
    { 
      bnpb[0][0]+=dbn;
      bnpb[1][1]-=.5*dbn;
      bnpb[2][2]-=.5*dbn;
    } 
    bnpb[0][1]=0; bnpb[1][0]=0; bnpb[0][2]=0; bnpb[2][0]=0;
    if (bug>0) printout("normal","bnpb %lg %lg %lg %lg %lg %lg\n",
                      bnpb[0][0],bnpb[1][1],bnpb[2][2],bnpb[0][1],bnpb[0][2],bnpb[1][2]); 
    
    Loop(j,0,3) Loop(k,0,3)
    { 
      bijw[j][k]=0;
      Loop(n,0,3) Loop(nn,0,3) bijw[j][k]+=wt[n][j]*wt[nn][k]*bnpb[n][nn];
    }
    
    if (bug>0) printout("normal","bijw %lg %lg %lg %lg %lg %lg\n",
                      bijw[0][0],bijw[1][1],bijw[2][2],bijw[0][1],bijw[0][2],bijw[1][2]); 
    Loop(n,0,3) bij[m+iall*n]=bijw[n][n];
    bij[m+iall*3]=bijw[0][1];
    bij[m+iall*4]=bijw[0][2];
    bij[m+iall*5]=bijw[2][1];
    
    if (whoelse[m]>0)  /* multiple point */
    { 
      Loop(k,1,whoelse[m][0]+1)
      { 
        j=whoelse[m][k];
        for (n=0;n<6;n++) bij[j+iall*n]=bij[m+iall*n];
      }
    }
  }
}
