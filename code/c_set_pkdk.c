/* contains c_set_pkdk */
#include "global.h"
/*  set P/k based on bij or a turbulent viscosity */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_set_pkdk(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep; /* needed arrays */
  char *cltp;
  double *q, *rho, *u[3], *bij=0, *visc=0;
  int *itrange; /* use if available */
  double *pkdk;  /* set */
  char *name; /* read 'bij' or name of turbulent viscosity */
  
  int i,j,k,its,ite,ip[4],i4dp[4],i4dm[4],ii[4],mc[8],ia[3],mmid,iallm,iall;
  double cave[8],rhom,qm,grad[8][3],volc,dudx[3][3],sij[3][3],ss,bm[3][3];
  int two[3]={2,2,2};
  int ij2n[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
  char borv='x';
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iall=Prod4(i4d);
  
  name=readname(fpin);
  printout("normal","set pkdk using %s",name);
  if (arraysize(name)==6*iall) 
  { borv='b'; bij=(double *)need(name); printout("normal"," for bij\n"); }
  else if (arraysize(name)==iallm)
  { borv='v'; visc=(double *)need(name); printout("normal"," for turb viscosity\n"); }
  else
  { 
    printout("error c_set_pkdk"," error %s is neither bij or turb viscosity, exit\n",name); 
    exitm4d(0); 
  }
  
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  q=(double *)need("qturb"); 
  rho=(double *)need("rho");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  
  pkdk=(double *)createarray("pkdk",iallm,'d',0);
  
  geomcinit();
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* each cont. c.v. */
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
    pkdk[mmid]=0;   /* initialize as zero */
    if (wherep[In4(ip,i4dp)]<0) continue;  /* not valid c.v. leave as zero */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    if (visc>0 && visc[mmid]==0) continue;
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }
    volc=geomcgradvol(ip,grad);
    Loop(i,0,3) Loop(j,0,3)   /*  velocity gradient */
    { 
      dudx[i][j]=0;
      Loop(k,0,8) dudx[i][j] += grad[k][j]*u[i][mc[k]];
    }
    fixdudx(dudx);    /* fix dudx for incompressible continuity */
    Loop(i,0,3) Loop(k,0,3) /* strain rate  */
    sij[i][k]=.5*(dudx[i][k]+dudx[k][i]);
    
    geomcvave(cave,ip); /* for average values */
    
    if (borv=='b') /* pkdk =-2bijSij */
    {
      qm=0; 
      Loop(i,0,3) Loop(j,0,3) bm[i][j]=0; 
      Loop(k,0,8)  
      { 
        qm+=cave[k]*q[mc[k]];
        Loop(i,0,3) Loop(j,0,3) bm[i][j] += cave[k]*bij[mc[k]+iall*ij2n[i][j]];
      }
      Loop(i,0,3) Loop(j,0,3) pkdk[mmid] += -2.*bm[i][j]*sij[i][j];
    }
    
    else if (borv=='v') /* pkdk=2SijSij*visc/(rho*q*q) */
    {
      rhom=0; qm=0;Loop(k,0,8)  
      { 
        qm+=cave[k]*q[mc[k]];
        rhom+=cave[k]*rho[mc[k]];
      }
      ss=0;  /* strain rate */
      Loop(i,0,3) Loop(k,0,3) ss += sij[i][k]*sij[i][k];
      ss*= 2;   
      pkdk[mmid]=ss*visc[mmid]/(qm*qm*rhom);
    }
  } /* each cont c.v. */
}	 
