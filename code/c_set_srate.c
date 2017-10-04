/* contains c_set_srate */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

/* set srate = strain rate  */

void c_set_srate(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep;   /* need */
  double *u[3]; 
  char *cltp;
  int *itrange; /* use if available */ 
  double *srate;  /* set */
  
  int i,j,k,i4dm[4],iallm,ip[4],i4dp[4],its,ite,ii[4],ia[3],mc[8],mmid,mwp;
  double grad[8][3],volc,dudx[3][3],sij[3][3],ss,s;
  int two[3]={2,2,2};
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  cltp=(char *)need("cltp");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1; }
  
  geomcinit();
  srate=(double *)createarray("srate",iallm,'d',0);
  
  /* loop over contiuity control volumes */
  Loop(ip[3],its,ite+1) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));  /* mid index */      
    srate[mmid]=0;  /* initialize as zero */
    mwp=wherep[In4(ip,i4dp)];
    if (mwp<0) continue;  /* not valid c.v. leave as zero */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    volc=geomcgradvol(ip,grad);
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }  
    Loop(i,0,3) Loop(j,0,3)   /*  velocity gradient */
    { 
      dudx[i][j]=0;
      Loop(k,0,8) dudx[i][j] += grad[k][j]*u[i][mc[k]];
    }
    fixdudx(dudx);    /* fix dudx for incompressible continuity */
    Loop(i,0,3) Loop(k,0,3) /* strain rate  */
    sij[i][k]=.5*(dudx[i][k]+dudx[k][i]);
    ss=0;  /* strain rate */
    Loop(i,0,3) Loop(k,0,3) ss += sij[i][k]*sij[i][k];
    ss*= 2;
    s=sqrt(ss);
    srate[mmid]=s;
  }
}
