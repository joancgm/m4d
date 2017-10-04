/* contains c_ddtall */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++) 
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* set an array of timescales ddtall[mmid,6] in the form 1/dt for between the points */

/*------------------ -------------------*/
void c_ddtall(FILE *fpin, FILE *fprint)
{  
  int *i4d,*wherep,*matchpc; /* need */
  double *U[3],*rho; 
  char *cltp;
  double *vlam=0,*zrotation=0;/* use if available */
  
  double *ddtall,*ddtmax; /* create */
  
  int i,j,k,i4dm[4],i4dp[4],iall,iallm,iallp,ip[4],mmid,ii[4],mc,ia[4];
  int two[3]={2,2,2};
  double volc,grad[8][3],ddt,dijk[3][3],Um[3],rhom;
  double rot[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double sij[3][3],wij[3][3],dudx[3][3],wijr[3][3],ss,wr,cmax,cc;
  
  /*
   n=0 convection (cfl) 
   n=1 laminar viscosity (using vlam if available)
   n=2 new totsafe
   n=3 vorticity/rotation parameter
   n=4 strain rate
   n=5 old MEFP+ rotation
   a companion turbulence ddt is om (from turb model, not set here)
   */
  int nc=0, nl=1, nt=2, nw=3, ns=4, nr=5, ntot=6;  /* corresponds to list above */
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  matchpc=(int *)need("matchpc");
  cltp=(char *)need("cltp");
  U[0]=(double *)need("U1");
  U[1]=(double *)need("U2");
  U[2]=(double *)need("U3");
  rho=(double *)need("rho");
  vlam=(double *)find("vlam"); 
  zrotation=(double *)find("zrotation");
  if (zrotation>0) 
  {  
    rot[0][1]=-zrotation[0]; 
    rot[1][0]=zrotation[0];
  }
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallp=Prod4(i4dp);
  iallm=Prod4(i4dm);
  
  ddtall=(double *)createarray("ddtall",iallm*ntot,'d',0);
  ddtmax=(double *)createarray("ddtmax",1,'d',0);
  Loop(i,0,iallm*ntot) ddtall[i]=0;
  ddtmax[0]=0;
  geomcinit();
  
  /* loop over contiuity control volumes */
  Loop(ip[3],0,i4d[3]) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    i=In4(ip,i4dp);
    if (wherep[i]<0) continue;   /* do only for non-zero, non-solid cont c.v. */
    if (cltp[i]=='s' || cltp[i]=='S') continue;
    
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3])); 	
    volc=geomcgradvol(ip,grad);   /* geometry arrays */
    
    Loop(j,0,3) Loop(k,0,3) { dijk[j][k]=0; dudx[j][k]=0; }
    Loop(j,0,3) Um[j]=0;
    rhom=0;
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      k=ia[0]+2*ia[1]+4*ia[2];
      mc=In4(ii,i4d);
      Loop(j,0,3) Um[j]+=.125*U[j][mc];
      rhom+=.125*rho[mc];
      Loop(i,0,3) Loop(j,0,3) dudx[i][j] += grad[k][j]*U[i][mc];  /*  velocity gradient */
    } 
    /* fixdudx(dudx); */   /* fix dudx for incompressible continuity */
    ss=0;
    wr=0;
    Loop(i,0,3) Loop(k,0,3) /* strain rate and vorticity */
    { 
      sij[i][k]=.5*(dudx[i][k]+dudx[k][i]);
      wijr[i][k]=.5*(dudx[i][k]-dudx[k][i]);
      wij[i][k]=wijr[i][k]+rot[i][k];
      ss+=sij[i][k]*sij[i][k];
      wr+=(rot[i][k]+wijr[i][k])*(rot[i][k]+wijr[i][k]);
    }
    ddtall[mmid+ns*iallm]=sqrt(2*ss);
    ddtall[mmid+nw*iallm]=sqrt(2*wr);
    Loop(j,0,3)
    {
      dijk[0][j]=.5*(-grad[0][j]+grad[1][j]-grad[2][j]+grad[3][j]
                     -grad[4][j]+grad[5][j]-grad[6][j]+grad[7][j]);
      dijk[1][j]=.5*(-grad[0][j]-grad[1][j]+grad[2][j]+grad[3][j]
                     -grad[4][j]-grad[5][j]+grad[6][j]+grad[7][j]);
      dijk[2][j]=.5*(-grad[0][j]-grad[1][j]-grad[2][j]-grad[3][j]
                     +grad[4][j]+grad[5][j]+grad[6][j]+grad[7][j]);
    }
    /* finish setting ddtall */
    ddt=0;
    Loop(i,0,3) 
    {
      ddt=max(ddt,abs(Dot(Um,dijk[i])));
      if (vlam>0) ddtall[mmid+nl*iallm]+=(vlam[mmid]/rhom)*Dot(dijk[i],dijk[i]);
    }
    ddtall[mmid+nc*iallm]=ddt; 
    ddtmax[0]=max(ddtmax[0],ddt);
    /* MEFP+ d/dt */
    cmax=0;
    Loop(i,0,3)
    {
      cc=0;
      /* MEFP + 
      cmax=max(cmax,abs(dudx[i][i]));
      Loop(k,0,3) cc+=(dudx[i][k]+2*rot[i][k])*(dudx[k][i]+2*rot[k][i]);
      cmax=max(cmax,sqrt(abs(cc))); */
      /* old MEFP  (with rot added ) */
       Loop(k,0,3) cc+=sqrt(abs((dudx[i][k]+2*rot[i][k])*(dudx[k][i]+2*rot[k][i])));
      cmax=max(cmax,cc);
    }
    ddtall[mmid+nr*iallm]=cmax;
    /* totsafe based on eqs 23,24 25 */
    cmax=0;
    Loop(i,0,3)
    {
      cc=0;
      Loop(k,0,3) cc+=(dudx[i][k]+2*rot[i][k])*(dudx[k][i]+2*rot[k][i]);
      if (cc>0)  cc=sqrt(abs(cc));
      else cc=0; 
      cmax=max(cc,cmax);
      /* cmax=max(cmax,abs(dudx[i][i])); */ /* eq 19 */
    }
    /* version 2 no sum - consider each separately */
  /*  cmax=0;
    Loop(i,0,3) Loop(k,0,3)
    {
      cc=(dudx[i][k]+2*rot[i][k])*(dudx[k][i]+2*rot[k][i]);
      if (cc>0) cmax=max(cmax,sqrt(abs(cc)));
    }
   */
    ddtall[mmid+nt*iallm]=cmax;
  }
  /* fill in approximation for zero cont volumes for plotting */
  Loop(ip[3],0,i4d[3]) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    i=In4(ip,i4dp);
    if (matchpc[i]<0) continue;
    iexpand(matchpc[i],i4dp,ia);
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3])); 
    j=ia[0]-1+i4dm[0]*(ia[1]-1+i4dm[1]*(ia[2]-1+i4dm[2]*ia[3]));
    if (j<0 || j>iallm-1) continue;
    Loop(k,0,ntot) ddtall[mmid+k*iallm]=ddtall[j+k*iallm];
  }
}
