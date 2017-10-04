/* contains c_coefinit */
#include "global.h"
/* initialize a generic coefficient set  with a selected number of coefficients
 Note: 0 is reserved for the set to be solved using eqnsolve.
 It's up to the user to keep track of the others by number */
#define Loop(n,a,b) for (n=a;n<b;n++)

void c_coefinit(FILE *fpin, FILE *fprint)
{
  int  *noindfwpts;    /* need */
  int *nocoefs, *coef_n, **coef_i; double **coef_c; /* create */
  int i,n;
  
  noindfwpts=(int *)need("noindfwpts");
  nocoefs=(int *)createarray("nocoefs",1,'i',0);
  nocoefs[0]=readint(fpin);
  n=noindfwpts[0];
  coef_n=(int *)createarray("coef_n",n,'i',0);
  coef_i=(int **)createarray("coef_i",n,'p',0);
  coef_c=(double **)createarray("coef_c",n,'p',0);
  Loop(i,0,n)
  { 
    coef_n[i]=0; coef_i[i]=0; coef_c[i]=0;
  }
  printout("normal","coef arrays created, number of coefs = %d\n",nocoefs[0]);
}


