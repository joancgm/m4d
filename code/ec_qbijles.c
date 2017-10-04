/* contains c_qbijles */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* determine parameters length dxb srate and bij to be used dynamically for LES 
 over a time step:   uiujLES = c L^2 S^2 dui/dt duj/dt /(2kles) where kles = duk/dt duk/dt
 dxb = sqrt(3/bddx) where bddx = 1/dx^2 + 1/dy^2 + 1/dz^2 (cartesian)
 set bij = dui/dt duj/dt /(2kles) -2kles/3 deltaij
 
 suggest but not set, L = min(.41y, .08dxb),  c=1.825^2  see les.more.notes.doc
  */

void ec_qbijles(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep;   /* need */
  double *u[3],*dudt[3]; 
  char *cltp;
  int *itrange; /* use if available */ 
  double *bijm, *dxb, *srate, *sijbij;  /* set */
  
  int i,j,k,i4dm[4],iallm,iall,ip[4],i4dp[4],its,ite,ii[4],ia[3],mc[8],mmid,mwp;
  double grad[8][3],volc,dudx[3][3],sij[3][3],ss,s,ddx[3][3];
  double dudxdt[3][3],dsijdt[3][3],dsdtdsdt,dsdt,dudtm[3],uiuj[6],twok;
  int two[3]={2,2,2};
  double bddx;
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iall=Prod4(i4d);
  
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  dudt[0]=(double *)need("dU1dt");
  dudt[1]=(double *)need("dU2dt");
  dudt[2]=(double *)need("dU3dt");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1; }
  
  geomcinit();
  bijm=(double *) createarray("bijm",6*iallm,'d',0);
  dxb=(double *)createarray("Lddx2",iallm,'d',0);
  srate=(double *)createarray("srate",iallm,'d',0);
  sijbij=(double *)createarray("sijbij",iallm,'d',0);
  
  /* loop over contiuity control volumes */
  Loop(ip[3],its,ite+1) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));  /* viscosity index */      
    dxb[mmid]=0;  /* initialize as zero */
    sijbij[mmid]=0;
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
    Loop(i,0,3)
    { dudtm[i]=0;
      Loop(k,0,8) dudtm[i]+=.125*dudt[i][mc[k]];
    }
    Loop(i,0,3) Loop(j,0,3)   /*  velocity gradient */
    {  
      dudx[i][j]=0;
      Loop(k,0,8) dudx[i][j] += grad[k][j]*u[i][mc[k]];
      dudxdt[i][j]=0;
      Loop(k,0,8) dudxdt[i][j] += grad[k][j]*dudt[i][mc[k]];
    }
    fixdudx(dudx);    /* fix dudx for incompressible continuity */
    Loop(i,0,3) Loop(k,0,3) /* strain rate  */
    sij[i][k]=.5*(dudx[i][k]+dudx[k][i]);
    ss=0;  /* strain rate */
    Loop(i,0,3) Loop(k,0,3) ss += sij[i][k]*sij[i][k];
    ss*= 2;
    s=sqrt(ss);
    srate[mmid]=s;
    Loop(i,0,3) Loop(k,0,3) /* strain rate  */
    dsijdt[i][k]=.5*(dudxdt[i][k]+dudxdt[k][i]);
    dsdtdsdt=0;  
    Loop(i,0,3) Loop(k,0,3) dsdtdsdt += dsijdt[i][k]*dsijdt[i][k];
    dsdtdsdt *= 2;
    dsdt=sqrt(dsdtdsdt);	
    Loop(j,0,3)   /* 1/dx each direction */
    { 
      ddx[0][j]=.5*(grad[1][j]-grad[0][j]+grad[3][j]-grad[2][j]
                    +grad[5][j]-grad[4][j]+grad[7][j]-grad[6][j]);
      ddx[1][j]=.5*(grad[2][j]-grad[0][j]+grad[3][j]-grad[1][j]
                    +grad[6][j]-grad[4][j]+grad[7][j]-grad[5][j]);
      ddx[2][j]=.5*(grad[4][j]-grad[0][j]+grad[5][j]-grad[1][j]
                    +grad[6][j]-grad[2][j]+grad[7][j]-grad[3][j]);			  
    }
    bddx=Dot(ddx[0],ddx[0])+Dot(ddx[1],ddx[1])+Dot(ddx[2],ddx[2])
    +2*(Dot(ddx[0],ddx[1])+Dot(ddx[0],ddx[2])+Dot(ddx[1],ddx[2]));
    if (bddx <=0)
    { 
      printout("warning ex_qbijles","error calculating dx, set to zero %lg\n",bddx);
      dxb[mmid]=0;
    }
    else dxb[mmid]=sqrt(3./bddx);
    Loop(k,0,6) uiuj[k]=0;
    Loop(k,0,8) 
    { 
      uiuj[0]+=.125*dudt[0][mc[k]]*dudt[0][mc[k]];
      uiuj[1]+=.125*dudt[1][mc[k]]*dudt[1][mc[k]];
      uiuj[2]+=.125*dudt[2][mc[k]]*dudt[2][mc[k]];
      uiuj[3]+=.125*dudt[0][mc[k]]*dudt[1][mc[k]];
      uiuj[4]+=.125*dudt[0][mc[k]]*dudt[2][mc[k]];
      uiuj[5]+=.125*dudt[1][mc[k]]*dudt[2][mc[k]];
    }
    twok=uiuj[0]+uiuj[1]+uiuj[2];
    for (i=0;i<6;i++) bijm[mmid+i*iallm]=0;
    if (twok>0) 
    { 
      /* alternate qmid[mmid]=cles*dxmin[mmid]*sqrts[dsdt]; */
      Loop(i,0,6) bijm[mmid+i*iallm]=uiuj[i]/twok;
      Loop(i,0,3) bijm[mmid+i*iallm] -=1./3.;
    }
    Loop(i,0,3) Loop(k,0,3)
    sijbij[mmid]=sij[0][0]*uiuj[0]+sij[1][1]*uiuj[1]+sij[2][2]*uiuj[2]
    +2*sij[0][1]*uiuj[3]+2*sij[0][2]*uiuj[4]+2*sij[1][2]*uiuj[5];
    sijbij[mmid]/=twok;
  }
}
