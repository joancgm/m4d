/* contains findex, find1dinterp, find2dinterp, interp3d, findex3d */
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Det(a,b,c) (*(c) *(*((a)+1) * *((b)+2) - *((a)+2) * *((b)+1)) +  *((c)+1) *( *((a)+2) * *(b) - *(a) * *((b)+2)) +   *((c)+2) *( *(a) * *((b)+1) - *((a)+1) * *(b) ))

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
/* ---------------------- */
int findex(double a, double *a1d, int id, double *f) 
/* range limited, returns initial value if all messed up */
{
  int iout,i; 
  *f=0; iout=0;
  if (a >=a1d[id-1]) {*f=1; iout=id-2;}
  else
  {
    for (i=0;i<id-1;i++)
      if (a<a1d[i+1])
      {
        iout=i; *f=(a-a1d[i])/(a1d[i+1]-a1d[i]);
        break;
      }
  }
  if (*f<0) *f=0;  if (*f>1) *f=1;
  /*  printout("normal",("findex %f    %f(1) %f(i) %f(i+1) %f(last) %d %f\n",
   a, array[0],array[iout],array[iout+1],array[id-1],i,f[0]);  */
  return iout;
}
/* ------------------------ */
void find1dinterp(double *x, double (*xp)[3], double *ff, int nylimit)
{ 
  double di[3],dj[3];
  int i;
  Loop(i,0,3)  {di[i]=x[i]-xp[0][i];   dj[i]=xp[1][i]-xp[0][i]; }
  ff[0]=Dot(di,dj)/Dot(dj,dj);
  if (nylimit>0)
  { if (ff[0]<0) ff[0]=0; else if (ff[0]>1) ff[0]=1; }
}
/* -------------------------*/
int find2dinterp(double *x, double (*xp)[2][3], double *ff, int nylimit, int itermax)
{ 
  double di[3],dj[3],dk[3],f[3],pt[3],df[3],dett,rhs[3],rhsm,rhsn,svi,svj,svk;
  int iter,n;
  f[0]=ff[0]; f[1]=ff[1]; f[2]=0;
  for (n=0;n<3;n++) 
  {  
    pt[n]=(1-f[0])*((1-f[1])*xp[0][0][n]+f[1]*xp[0][1][n])
    +f[0]*((1-f[1])*xp[1][0][n]+f[1]*xp[1][1][n]);
    rhs[n]=x[n]-pt[n];
  }
  for (iter=0; iter<itermax; iter++)
  { 
    for (n=0;n<3;n++)
    { 
      di[n]=(1-f[1])*(xp[1][0][n]-xp[0][0][n])+f[1]*(xp[1][1][n]-xp[0][1][n]);
      dj[n]=(1-f[0])*(xp[0][1][n]-xp[0][0][n])+f[0]*(xp[1][1][n]-xp[1][0][n]);
    }
    Cross(di,dj,dk);
    dett=Dot(dk,dk);
    rhsm=Dot(rhs,rhs);
    if (dett==0) return(-1);
    svi=di[0]; svj=dj[0];svk=dk[0];
    di[0]=rhs[0]; dj[0]=rhs[1]; dk[0]=rhs[2];
    df[0]=Det(di,dj,dk);
    di[0]=svi; dj[0]=svj; dk[0]=svk;
    di[1]=rhs[0]; dj[1]=rhs[1]; dk[1]=rhs[2];
    df[1]=Det(di,dj,dk);
    f[0]+=df[0]; f[1]+=df[1];
    if (nylimit>0) 
    { 
      if (f[0]<0) f[0]=0; 
      if (f[0]>1) f[0]=1; 
      if(f[1]<0) f[1]=0; 
      if(f[1]>1) f[1]=1;
    }
    for (n=0;n<3;n++) 
    {  
      pt[n]=(1-f[0])*((1-f[1])*xp[0][0][n]+f[1]*xp[0][1][n])
      +f[0]*((1-f[1])*xp[1][0][n]+f[1]*xp[1][1][n]);
      rhs[n]=x[n]-pt[n];
    }
    rhsn=Dot(rhs,rhs);
    if (rhsn>=rhsm) return (iter);
    ff[0]=f[0]; ff[1]=f[1];
  }
  return(iter);
}
/* -------------------------- */
void interp3d(double *x, int *idim, double *f, double *xout) 
/* f is the non-interger index starting at 0, x dims must be >=2*/
{
  int i,ii[3];
  double ff[3];
  for (i=0;i<3;i++) 
  {
    ii[i]=f[i]; if (ii[i]>idim[i]-2) ii[i]=idim[i]-2;
	 ff[i]=f[i]-ii[i];
  }
  /*printout("normal",("interp3d idim %d %d %d, f %lg %lg %lg d ff %lg %lg %lg ii %d %d %d ",
   idim[0],idim[1],idim[2],f[0],f[1],f[2],ff[0],ff[1],ff[2],ii[0],ii[1],ii[2]); */
  *xout=(1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] 
  +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] 
  +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] 
  +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] 
  + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] 
  +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] 
  +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] 
  +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))];
}
/*-------------------------------------*/
/* find 3d interpolation index of xx in array xyz dimensioned 'id'*3   
 limit index range from ilo to ihi  tol is in terms of ff */
int findex3d(double *xx, double **x, int *id, int *ilo, int *ihi, 
             double tol, double *fout)
{
  int i,itmax,iter,ib[3][3],ia[3],nych=0,iamax[3];
  double f[3],xnow[3],distn,t,df[3],ft[3],distt;
  Loop(i,0,3) /* initial limits*/
  { 
    f[i]=fout[i];
	 if (ihi[i]>id[i]-1) ihi[i]=id[i]-1;
	 if (f[i]<ilo[i] || f[i]>ihi[i]) { f[i]=ilo[i]+ihi[i]; f[i]*=.5; }
	 df[i]=.3*(ihi[i]-ilo[i]);
	 iamax[i]=3; if (ilo[i]==ihi[i]) iamax[i]=1;
  }
  Loop(i,0,3) interp3d(x[i],id,f,&xnow[i]);   /* current point and distance */
  distn=(xnow[0]-xx[0])*(xnow[0]-xx[0])+(xnow[1]-xx[1])*(xnow[1]-xx[1])+
  (xnow[2]-xx[2])*(xnow[2]-xx[2]);
  /*printout("normal",("findex3d f %lg %lg %lg df %lg %lg %lg distn %lg\n",f[0],f[1],f[2],df[0],df[1],df[2],distn); */
  if (distn==0) return(0);
  
  t=1; Loop(i,0,3) if (t<ihi[i]-ilo[i]+1) t=ihi[i]-ilo[i]+1;
  itmax=10*(1+log(t/tol)/log(2.));
  Loop(i,0,3) {ib[i][0]=0; ib[i][1]=-1; ib[i][2]=1; }
  
  Loop(iter,0,itmax)
  {
    Loop3(ia,0,iamax)   /* find better point */
    {
      nych=0;
      if (ia[0]+ia[1]+ia[2]==0) continue;
      Loop(i,0,3) 
      {
        ft[i]=f[i]+df[i]*ib[i][ia[i]];
        if (ft[i]<ilo[i]) ft[i]=ilo[i];
        if (ft[i]>ihi[i]) ft[i]=ihi[i];
      } 
      Loop(i,0,3)  interp3d(x[i],id,ft,&xnow[i]);    /* test  point and distance */
      distt=(xnow[0]-xx[0])*(xnow[0]-xx[0])+(xnow[1]-xx[1])*(xnow[1]-xx[1])+
      (xnow[2]-xx[2])*(xnow[2]-xx[2]);
      /*	 printout("normal",("ft %lg %lg %lg  distt %lg distn %lg\n",ft[0],ft[1],ft[2],distt,distn); */
      if (distt<distn)  /* better point */
      { 
        Loop(i,0,3) f[i]=ft[i];
        distn=distt;
        Loop(i,0,3)
        {
          if (ia[i]==2) {ib[i][1]*=-1; ib[i][2]*=-1;}
          else if (ia[i]==1)  df[i]*=1.3;
        }
        nych=1;
        break;
      }
    }
	 if (nych==0) 
    { 
      if (df[0]<=tol && df[1]<=tol && df[2]<=tol) break;
      Loop(i,0,3) df[i]*=.5;
    }
  }
  Loop(i,0,3) fout[i]=f[i];
  
  return (iter); 
}


