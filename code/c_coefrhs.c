/* contains c_coefrhs coefrhs */
#include "global.h"
/* calculate the rhs of an equation based on the coefs at jcoef */
/* if roundoffd exists evaluate with roundoff filter */

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_coefrhs(FILE *fpin, FILE *fprint)
{ 
  char *namev,*namerhs;
  int jcoef;
  double fac; char oldnew;
  namev=readname(fpin);
  namerhs=readname(fpin);
  oldnew=read1charname(fpin);
  jcoef=readint(fpin);
  fac=readdouble(fpin);
  printout("normal","rhs for %s, %s, updated(-=) with %lg*coefs %d, oldnew=%c\n",
          namev,namerhs,fac,jcoef,oldnew);
  coefrhs(fprint,namev,namerhs,jcoef,oldnew,fac);
}
/* -------------------------- */
void coefrhs(FILE *fprint, const char *namev, const char *namerhs, 
             int jcoef, char oldnew, double fac)
{
  int *noindfwpts, *nogpts, *coef_n, **coef_i, *nocoefs;   /* needed arays */
  int *itrange; double *roundoff;/* use if available */
  double **coef_c, *v,round;
  double *rhs,drhs; /* set */
  
  int i,j,k,ns,ne,irep,neqs,npts;
  
  nocoefs=(int *)need("nocoefs");
  if (jcoef>nocoefs[0]-1)
  { 
    printout("error coefrhs"," error, coef only dimensioned for %d coefs\n",nocoefs[0]);
    exitm4d(0);
  }
  if (oldnew!='n' && oldnew!='o') 
  { 
    printout("error coefrhs","error, oldnew %c should be o or n\n",oldnew); 
    exitm4d(0); 
  }
  noindfwpts=(int *)need("noindfwpts");
  neqs=noindfwpts[0];
  nogpts=(int *)need("nogpts");
  npts=nogpts[0];
  itrange=(int *)find("itrange");
  roundoff=(double *)find("roundoffd");
  if (itrange==0) {ns=0,ne=neqs;}
  else {ns=noindfwpts[itrange[0]+1]; ne=noindfwpts[itrange[1]+2]; }
  coef_n=(int *)need("coef_n");
  coef_i=(int **)need("coef_i");
  coef_c=(double **)need("coef_c");
  v=(double *)need(namev);
  irep=arraysize(namev)/npts;
  if (irep>1) printout("normal","irep %d\n",irep); 
  rhs=(double *)find(namerhs);
  if (rhs==0)
  {	
    rhs=(double *)createarray(namerhs,neqs*irep,'d',0);
    Loop(i,0,neqs*irep) rhs[i]=0;
    printout("normal"," %s created\n",namerhs);
  }
  else if (oldnew=='n') Loop(i,ns,ne) Loop(k,0,irep) rhs[i+k*neqs]=0;
  
  Loop(i,ns,ne)
  if (coef_n[i]>0)
  { 
    Loop(k,0,irep)
    {
      drhs=0;
    Loop(j,0,coef_n[i]) 
      drhs-=fac*coef_c[i][j+coef_n[i]*jcoef]*v[coef_i[i][j]+k*npts];
    if (roundoff>0 && roundoff[0]>0)
    {
      round=0;
      Loop(j,0,coef_n[i]) 
      round+=abs(fac*coef_c[i][j+coef_n[i]*jcoef]*v[coef_i[i][j]+k*npts]);
      if (abs(drhs)<roundoff[0]*round) drhs=0;
    }
      rhs[i+k*neqs]+=drhs;
    }
  }
}


