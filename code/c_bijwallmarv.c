/* contains c_bijwallmarv */
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Loop(n,a,b) for (n=a;n<b;n++)

void c_bijwallmarv(FILE *fpin, FILE *fprint)  

/* set bij at walls consistent with MA/MARV model */
{
  int   *i4d,*noindwpts,*whoisw,*wnear,**whoelse;     /* need */
  double  *wnorm,*u[3],*om,*xyz; 
  int *itrange; /* use if available */
  double *bij; /* modify */ 
  
  int its,ite;
  int iall,i,j,k,m,n,inw;
  double up[3],usub,*wn,umag,dx[3],dy,ub[3],fexp,bnpb[3],bijw[6];
  
  i4d=(int *)need("idim4d");             /* find input values and arrays */
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  itrange=(int *)find("itrange");
  if (itrange>0) {its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1;}
  wnear=(int *)need("wnear");
  noindwpts=(int *)need("noindwpts");
  whoisw=(int *)need("whoisw");
  whoelse=(int **)need("whoelse");
  xyz=(double *)need("xyz");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  om=(double *)need("omturb");
  wnorm=(double *)need("wnorm");
  bij=(double *)need("bij");
  
  Loop(i,noindwpts[its+1],noindwpts[ite+2])
  { 
    m=whoisw[i];
    inw=wnear[i];
    wn=wnorm+3*i;
    Loop(n,0,3) up[n]=u[n][inw];
    usub=Dot(up,wn);
    Loop(n,0,3) up[n]-=usub*wn[n];
    umag=sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
    Loop(n,0,3) dx[n]=xyz[inw+n*iall]-xyz[m+n*iall];
    dy=abs(Dot(wn,dx));
    if (umag>.001*abs(usub) && dy>0 && om[m]>0)   /* ok to do this */
    { 
      Loop(n,0,3) up[n]=up[n]/umag;  /* flow paralled vector */
      Cross(wn,up,ub); /* binormal vector */
      bnpb[0]=-.33333;   /* normal */
      fexp=.0384*umag/(dy*om[i]);
      if (fexp>.0001) fexp=1.-exp(-fexp);
      bnpb[2]=1./6.-(.2+1./6.)*fexp; /* binormal */
      bnpb[1]=-bnpb[0]-bnpb[2];  /* parallel */
    }
    else  /* insufficient information to set bij wall, assume zero velocity*/
    { 
      up[0]=1; up[1]=0; up[2]=0;
      usub=Dot(wn,up);
      if (abs(usub)>.999) {up[0]=0; up[1]=1; usub=Dot(wn,up);}
      Loop(n,0,3) up[n]-=usub*wn[n];
      umag=sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
      Loop(n,0,3) up[n]/=umag;
      Cross(wn,up,ub);
      bnpb[0]=-.33333;
      bnpb[1]=-.5*bnpb[0];
      bnpb[2]=bnpb[1];
    }
    bijw[0]=bnpb[0]*wn[0]*wn[0]+bnpb[1]*up[0]*up[0]+bnpb[2]*ub[0]*ub[0];
    bijw[1]=bnpb[0]*wn[1]*wn[1]+bnpb[1]*up[1]*up[1]+bnpb[2]*ub[1]*ub[1];
    bijw[2]=bnpb[0]*wn[2]*wn[2]+bnpb[1]*up[2]*up[2]+bnpb[2]*ub[2]*ub[2];
    bijw[3]=bnpb[0]*wn[0]*wn[1]+bnpb[1]*up[0]*up[1]+bnpb[2]*ub[0]*ub[1];
    bijw[4]=bnpb[0]*wn[0]*wn[2]+bnpb[1]*up[0]*up[2]+bnpb[2]*ub[0]*ub[2];
    bijw[5]=bnpb[0]*wn[1]*wn[2]+bnpb[1]*up[1]*up[2]+bnpb[2]*ub[1]*ub[2];
    for (n=0;n<6;n++) bij[m+iall*n]=bijw[n];
    if (whoelse[m]>0)  /* multiple point */
    { 
      Loop(k,1,whoelse[m][0]+1)
      { 
        j=whoelse[m][k];
        for (n=0;n<6;n++) bij[j+iall*n]=bijw[n];
      }
    }
  }
}
