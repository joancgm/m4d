/* contains omwall, c_omwallcoakley, c_omwallmarv, c_omwallmarvs, c_omwall */
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Loop(n,a,b) for (n=a;n<b;n++)

/* --------------------------- */
/* set om at walls = com dq/dy */
void omwall(FILE *fprint, double com)  
{
  int   *i4d,*noindwpts,*whoisw,*wnear,**whoelse;     /* need */
  double  *wnorm,*q,*xyz; 
  int *itrange; /* use if available */
  double *om; /* modify */ 
  
  int its,ite;
  int iall,i,k,m,n,inw;
  double *wn,dx[3],dy;
  
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
  om=(double *)need("omturb");
  wnorm=(double *)need("wnorm");
  q=(double *)need("qturb");
  
  Loop(i,noindwpts[its+1],noindwpts[ite+2])
  {
    m=whoisw[i];
    inw=wnear[i];
    wn=wnorm+3*i;
    Loop(n,0,3) dx[n]=xyz[inw+n*iall]-xyz[m+n*iall];
    dy=abs(Dot(wn,dx));
    om[m]=com*q[inw]/dy;
    if (whoelse[m]>0)  /* multiple point */
      Loop(k,1,whoelse[m][0]+1) om[whoelse[m][k]]=om[m];
  }
}
/* -------------------------------- */
/* set om at walls consistent with modified Coakley model (vki notes) */
void c_omwallcoakley(FILE *fpin, FILE *fprint)  
{ 
  omwall(fprint,.56);
}
/* --------------------------------- */
/* set om at walls consistent with MA/MARV model */
void c_omwallmarv(FILE *fpin, FILE *fprint)  
{   
  omwall(fprint,.089);
}/* --------------------------------- */
/* set om at walls consistent with MARVS model */
void c_omwallmarvs(FILE *fpin, FILE *fprint)  
{   
  omwall(fprint,.095);
}
/* --------------------------------- */               
/* set om at walls  om(b+c)= b a dq/dy + c omnw  see note omwall.doc!!   read a,b,c */
void c_omwall(FILE *fpin, FILE *fprint)  
{
  int   *i4d,*noindwpts,*whoisw,*wnear,**whoelse;     /* need */
  double  *wnorm,*q,*xyz,*vlamg; 
  int *itrange; /* use if available */
  double *om; /* modify */ 
  double a,b,c,cc;
  
  int its,ite;
  int iall,i,k,m,n,inw;
  double *wn,dx[3],dy;
  
  a=readdouble(fpin);
  b=readdouble(fpin);
  c=readdouble(fpin);
  
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
  om=(double *)need("omturb");
  wnorm=(double *)need("wnorm");
  q=(double *)need("qturb");
  vlamg=(double *)need("vlamg");
  
  Loop(i,noindwpts[its+1],noindwpts[ite+2])
  {
    m=whoisw[i];
    inw=wnear[i];
    wn=wnorm+3*i;
    Loop(n,0,3) dx[n]=xyz[inw+n*iall]-xyz[m+n*iall];
    dy=abs(Dot(wn,dx));
    cc=2*c*sqrt(vlamg[m])/dy;
    om[m]=(b*sqrt(a*q[inw]/dy)+cc)/(b+cc/sqrt(om[inw]));
    om[m]=om[m]*om[m];		 
    if (whoelse[m]>0)  /* multiple point */
      Loop(k,1,whoelse[m][0]+1) om[whoelse[m][k]]=om[m];
  }
}
