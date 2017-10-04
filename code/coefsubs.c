/* contains coefcadd, coefcombine, coefcenterfirst, coeforder, coeforderv */
#include "global.h" 
#define Loop(n,a,b) for (n=a;n<b;n++)

/* --------------coefficient service subroutines -------------*/

/* add n points ipt with coef values cadd to the jcoef set of coef's in cc */
void coefcadd(int iw, int *cn, int **ci, double **cc, int nocoefs, int jcoef, int n, int *ipt, double *cadd)
{ 
  int i,j,nnew,*cinew,*newai,n1,nyok;
  double *ccnew;
  /* printout("normal","coefcadd iw %d nocoefs %d jcoef %d n %d cn %d\n",iw,nocoefs,jcoef,n,cn[iw]);
	Loop(i,0,n) printout("normal"," %d %lg,",ipt[i],cadd[i]); printout("normal","\n"); */
  if (iw<0) 
  {
    printout("warning coefcadd","\n *** warning coefcadd call with neg eq No. ignored\n");
    return;
  }
  newai=(int*)smalloca(n,'i');
  nnew=0;
  n1=cn[iw];
  if (n1==0) 
	 Loop(i,0,n) {if (cadd[i]!=0) {newai[nnew]=i; nnew++; }}
  else
	 Loop(i,0,n) 
  {
    if (cadd[i]!=0)
    {
      nyok=0; 
      Loop(j,0,n1) if (ci[iw][j]==ipt[i]) 
      {
        nyok=1; 
        cc[iw][j+n1*jcoef]+=cadd[i]; 
        break; 
      }
      if (nyok==0) {newai[nnew]=i; nnew++; }
    }
  }
  if (nnew>0)
  {
    cinew=(int *)smalloca(cn[iw]+nnew,'i');
    ccnew=(double *)smalloca((cn[iw]+nnew)*nocoefs,'d');
    if (n1>0) 
    {
      Loop(i,0,cn[iw])
      {
        cinew[i]=ci[iw][i];
        Loop(j,0,nocoefs) ccnew[i+j*(nnew+n1)]=cc[iw][i+n1*j];
      }
      free(ci[iw]);  free(cc[iw]); 
    }
    ci[iw]=cinew; cc[iw]=ccnew;
    cn[iw]+=nnew;
    Loop(i,0,nnew)
    {
      ci[iw][n1+i]=ipt[newai[i]]; 
      Loop(j,0,nocoefs) cc[iw][n1+i+j*(nnew+n1)]=0; 
      cc[iw][n1+i+jcoef*(nnew+n1)]=cadd[newai[i]];
    }
  }
  free(newai);
}

/*---------------------------------------------------*/
/* combine coeffcients with the same indices */
void coefcombine(int iw, int *cn, int **ci, double **cc, int irep, double tol)
{
  int nn,n,i,j,k,*newai; double *cnew,cmax[irep],cm;
  nn=cn[iw]; 
  if (nn<2) return;
  k=nn;
  Loop(j,0,irep) cmax[j]=abs(cc[iw][j*nn]);
  Loop(n,1,nn)    /* combine coef */
  { 
    Loop(j,0,irep) cmax[j]=max(cmax[j],abs(cc[iw][n+j*nn]));
    Loop(i,0,n) 
    {
      if (ci[iw][n]==ci[iw][i])
      {
        Loop(j,0,irep) cc[iw][i+j*nn]+=cc[iw][n+j*nn];
        ci[iw][n]=-1;
        k--; break;
      }
    }
  }
  Loop(n,0,nn) /* clean off roundoff */
  Loop(j,0,irep) {if (abs(cc[iw][n+j*nn])<tol*cmax[j]) cc[iw][n+j*nn]=0; }
  Loop(n,0,nn)  if (ci[iw][n]!=-1) /* clean off zero points */
  { 
    cm=0;
	 Loop(j,0,irep) cm=max(cm,abs(cc[iw][n+j*nn]));
	 if (cm==0) { ci[iw][n]=-1; k--;}
  }
  
  if (k==0) {cn[iw]=0;  free(ci[iw]); free(cc[iw]); ci[iw]=0; cc[iw]=0; }
  
  else if (k<nn)
  {
    newai=(int *) smalloca(k,'i');
	 cnew=(double *)smalloca(k*irep,'d');
	 i=0;
	 Loop(n,0,nn) if (ci[iw][n]>=0)
    { 
      newai[i]=ci[iw][n];
      Loop(j,0,irep) cnew[i+j*k]=cc[iw][n+j*nn]; 
      i++;
    } 
    free(ci[iw]);  ci[iw]=newai; 
    free(cc[iw]);    cc[iw]=cnew; 
    cn[iw]=k;  
  }  
}
/*-------------------------------------------*/
/* if there are any coefficients, make the center one first */
void coefcenterfirst(int icen, int iw, int *cn, int **ci, double **cc, int jrep)
{ 
  int i,j,k,nn,*ii,is; double *c,cs;
  nn=cn[iw];  
  if (nn==0) return;
  ii=ci[iw]; c=cc[iw]; 
  Loop(i,0,nn)  { if (ii[i]==icen) break;}
  if (i<nn)
  {
    if (i!=0) 
    { 
      is=ii[0]; ii[0]=ii[i]; ii[i]=is;
      Loop(j,0,jrep)
      {
        cs=c[0+j*nn]; c[0+j*nn]=c[i+j*nn]; c[i+j*nn]=cs;
      }
    }
  }
  else  /* add center point coefficient */
  {
    ii=(int *)smalloca((nn+1),'i');
	 c=(double *)smalloca((nn+1)*jrep,'d');
	 Loop(k,0,nn)
    { 
      ii[k+1]=ci[iw][k];
      Loop(j,0,jrep) c[k+1+j*(nn+1)]=cc[iw][k+j*nn];
    }
	 ii[0]=icen;
	 Loop(j,0,jrep) c[j*(nn+1)]=0;
	 free(ci[iw]); free(cc[iw]);
	 ci[iw]=ii; cc[iw]=c; 
	 nn++; cn[iw]=nn;
  }
}
/*---------------------------------------------*/
/* call centerfirst then order remianing coefficients  in increasing index */

void coeforder(int icen, int iw, int *cn, int **ci, double **cc, int jrep)
{
  int i,j,k,nn,*ii,is,kmin,kat; double *c,cs;
  
  coefcenterfirst(icen,iw,cn,ci,cc,jrep);
  nn=cn[iw];  
  if (nn==0) return;
  ii=ci[iw]; c=cc[iw]; 
  
  Loop(k,1,nn)  /* set k to lowest coef left */
  { 
    kat=k; kmin=ii[k];
	 Loop(i,(k+1),nn)
    if (ii[i]<kmin) {kmin=ii[i]; kat=i;}
	 if (kat!=k)
    {
      is=ii[k]; ii[k]=ii[kat]; ii[kat]=is;
      Loop(j,0,jrep)
      {cs=c[k+j*nn]; c[k+j*nn]=c[kat+j*nn]; c[kat+j*nn]=cs; }
    }
  }
}
/*--------------------------------------*/
/* call centerfirst then order remianing coefficients  in decreasing value */

void coeforderv(int icen, int iw, int *cn, int **ci, double **cc, int jrep, int jorder)
{ 
  int i,j,k,nn,*ii,is,kmin,kat; double *c,cs,cmax;
  
  coefcenterfirst(icen,iw,cn,ci,cc,jrep);
  nn=cn[iw];  
  if (nn==0) return;
  ii=ci[iw]; c=cc[iw]; 
  
  Loop(k,1,nn)  /* set k to  highest coef left*/
  {
    kat=k; kmin=ii[k]; cmax=c[k+jorder*nn];
	 Loop(i,(k+1),nn)
    if (c[i+jorder*nn]>cmax) {cmax=c[i+jorder*nn]; kat=i;}
	 if (kat!=k)
    { 
      is=ii[k]; ii[k]=ii[kat]; ii[kat]=is;
      Loop(j,0,jrep)
      {cs=c[k+j*nn]; c[k+j*nn]=c[kat+j*nn]; c[kat+j*nn]=cs; }
    }
  }
}



