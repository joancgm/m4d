/* contains  eqnrhs, eqnrhsbij, eqncplus, eqnerror, eqnerrorwho, eqnbndry, eqncenter, eqncenterbij */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
/* ------------------------------ */
void eqnrhs(int ns, int nn, int *cn, int **ci, double **cc, double *ch, 
            double *rhsin, double *rhs)
{
  int ieq,i;
  
  Loop(ieq,ns,nn)
  { 
    rhs[ieq]=rhsin[ieq];
    if (cn[ieq]>0)
      Loop(i,0,cn[ieq])
      rhs[ieq]-=ch[ci[ieq][i]]*cc[ieq][i];
  }
}
/* ----------------------------- */
void eqnrhsbij(int ns, int nn, int *cn, int **ci, double **cc, int npts, 
               int neqs, double *beqm, double *dbij, double *beqd,
               char *nybijeq, double *rhs)
{
  int ieq,i,j,k; int n=6; /* if more general routine is needed make n a parameter */
  /* note first value in ci must be the point itself */
  int nobug=0,	buga[2]={1,10}, bugb[2]={11,162};  /* eq range for print */
  Loop(ieq,ns,nn)
  {
    Loop(j,0,n) rhs[ieq+j*neqs]=0;
    if (cn[ieq]>0)
    {
      Loop(j,0,n) 
      if (nybijeq[ieq+neqs*max(j-2,0)]=='y')
      {
        rhs[ieq+j*neqs]=beqm[ieq+j*neqs];			 
        Loop(i,0,cn[ieq])
        {
          rhs[ieq+j*neqs]-=cc[ieq][i]*dbij[ci[ieq][i]+j*npts];
          if (nobug>0)
            Loop(k,0,nobug)
            if (ieq>buga[k] && ieq<bugb[k])
              printout("normal","ieq %d j %d beqm %lg dbij %d %lg cc %lg rhs %lg\n",
                     ieq,j,beqm[ieq+j*neqs],ci[ieq][i],
                     dbij[ci[ieq][i]+j*npts],cc[ieq][i],
                     rhs[ieq+j*neqs]);
        }
        Loop(i,0,n)  
        {
          rhs[ieq+j*neqs]+=beqd[ieq+neqs*(j+n*i)]*(dbij[ci[ieq][0]+i*npts]);
        }
      }
    }
  }
}

/*----------------------*/
void eqncplus(int ns, int nn, int *cn, double **cc, double *cplus)
{ 
  int ieq,i;
  
  Loop(ieq,ns,nn)
  { 
    cplus[ieq]=0;
    if (cn[ieq]>0)
      Loop(i,0,cn[ieq])
      if (cc[ieq][i]>0) cplus[ieq]+=cc[ieq][i];
  }
}
/*----------------------*/
void eqnerror(int ns, int nn, double *rhs, double *cplus, double *err, int ierrtot)
{ 
  int i,j,nactive;
  double emin=0,emax=0,esq,r;
  esq=0; nactive=0;
  Loop(i,ns,nn)
  {
    if (cplus[i]>0)
    { 
      r=rhs[i]/cplus[i]; 
		if (nactive==0) {emin=r; emax=r; }
		else {emin=min(emin,r); emax=max(emax,r); }
		esq+=r*r;
		nactive++;
    }
  }
  esq=sqrt(esq/nactive);
  for (i=ierrtot-1;i>0;i--)    /* move down for more values */
    Loop(j,0,3) err[j+3*i]=err[j+3*(i-1)];
  err[0]=esq; err[1]=emin; err[2]=emax;
}
/*---------return nactive; in err: rms,min,max; in im,ix; eq nos for min and max */
int eqnerrorwho(int ns, int nn, double *rhs, double *cplus, char *ny,
                double *err, int *im, int *ix)
{ 
  int i,nactive=0,imin=0,imax=0;
  double emin=0,emax=0,esq=0,r;
  Loop(i,ns,nn)
  if (cplus[i]>0 &&  (ny==0 || ny[i]=='y') )
  { 
    r=rhs[i]/cplus[i];  
    if (nactive==0) {emin=r; emax=r; imin=i; imax=i;}
    else if (r<emin) {emin=r; imin=i;}
    else if (r>emax) {emax=r; imax=i;}
    esq+=r*r;
    nactive++;
  }
  if (nactive>0) esq=sqrt(esq/nactive);
  err[0]=esq; err[1]=emin; err[2]=emax;
  *im=imin; *ix=imax;
  return nactive;
}
/*----------------------------*/
void eqnbndry(int *ipeqexit, int *cn, int **ci, double **cc, double *rhsin, double *rhs, double *ch)
{
  int ii,i; 
  if (ipeqexit==0) return;
  Loop(ii,0,ipeqexit[0]) 
  { i=ipeqexit[ii+1];
    if (cn[i]>0) 
    {
      if (cc[i][0]>0) 
      { 
        eqnrhs(i,i+1,cn,ci,cc,ch,rhsin,  rhs);
        ch[ci[i][0]]+=rhs[i]/cc[i][0]; 
      }
    }
  }
}
/*----------------------------*/
void eqncenter(int ns, int nn, double *cplus, double *rhs, int **ci, double *ch)
{ 
  int i;
  Loop(i,ns,nn) if (cplus[i]>0) ch[ci[i][0]]+=rhs[i]/cplus[i]; 
}

/*---------------------           -------*/
void eqncenterbij(int ns, int nn, int neqs, int npts, double *cplus, double *beqd, double *rhs, int **ci, char *nybijeq, double *bij, double *dbij)
/* err anal returned in ddberr: max ddb, rms ddb, no pts with real reduction */
{ 
  int i,j,k,n,ia,ja,ka;
  double aa[6][6],bb[6],out[6],ah,bs,as;
  int	buga=-1, bugb=-1;   /* eq no range to print coef and results*/
  
  Loop(i,ns,nn) if (nybijeq[i]=='y') 
  { 
    Loop(j,0,6) Loop(k,0,6) aa[j][k]=0;  /* set up 6x6 coefficients */
    aa[0][0]=1; aa[0][1]=1; aa[0][2]=1; bb[0]=0;  /* 11 equation from -22-33 */
    Loop(j,0,3) bb[0]-=bij[ci[i][0]+npts*j]+dbij[ci[i][0]+npts*j];
    Loop(j,1,6) 
    {
      if (nybijeq[i+neqs*max(j-2,0)]=='y')
		{ 
        aa[j][j]+=cplus[i];
        bb[j]=rhs[i+neqs*j];
        Loop(k,0,6) aa[j][k]-=beqd[i+neqs*(j+6*k)];
		}
		else {aa[j][j]=1; bb[j]=0; }
    }
    if (i>buga && i< bugb) Loop(j,0,6)
    {
      printout("normal","eq %d j %d bb %lg cplus %lg\n",i,j,bb[j],cplus[i]);
      Loop(k,0,6) printout("normal"," %lg",aa[j][k]);
      printout("normal","\n");
    }
    
    for (ia=0;ia<6;ia++) /* solve */
    {
      if (ia>0)
		{ 
        ah=abs(aa[ia][ia]); ja=ia;  /* find best equation */
        for (n=ia+1;n<6;n++) if (abs(aa[n][ia])>ah) {ah=abs(aa[n][ia]); ja=n;}
        if (ja!=ia)             /* switch equations */
        {
          bs=bb[ia]; bb[ia]=bb[ja]; bb[ja]=bs;
          for (ka=0;ka<n;ka++) 
          { as=aa[ia][ka]; aa[ia][ka]=aa[ja][ka]; aa[ja][ka]=as; }
        } /* ja!=ia */
		}
      bs=aa[ia][ia];     /* normalize */
      bb[ia] /= bs;
      for (n=ia;n<6;n++) aa[ia][n] /= bs;
      for (n=ia+1;n<6;n++)   /* reduce */
      { 
        bs= aa[n][ia]; bb[n] -= bs*bb[ia];
        for (ja=ia;ja<6;ja++) aa[n][ja] -= bs*aa[ia][ja];
      }	    
    } /* ia */
    out[5]=bb[5];   /* back substitute */
    for (ia=4;ia>=0;ia--)
    { 
      out[ia]=bb[ia];
      for (ja=ia+1;ja<6;ja++) out[ia] -= out[ja]*aa[ia][ja]; 
    } 
    if (i>buga && i< bugb) 
    {printout("normal","out"); Loop(k,0,6) printout("normal"," %lg",out[k]); printout("normal","\n"); }
    
    Loop(j,0,6) dbij[ci[i][0]+j*npts]+=out[j];
  }
}
