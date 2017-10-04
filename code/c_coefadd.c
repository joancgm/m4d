/* contains c_coefadd */
#include "global.h"
/* form a sum of coef sets  */

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_coefadd(FILE *fpin, FILE *fprint)
{ 
  int *noindfwpts,*cn,*nocoefs;/* needed arrays */
  double **cc;  /* modify */
  int *itrange; /* use if available */
  int jto, nfrom, *jfrom;  double *fac; /* input parameters */
  
  int i,ieq,j,neqs,neqe,jmax;
  double sum;
  
  jto=readint(fpin);
  jmax=jto;
  nfrom=readint(fpin);
  jfrom=(int *)tmalloca(nfrom,'i');
  fac=(double *)tmalloca(nfrom,'d');
  printout("normal","coef_c.%d=",jto);
  Loop(i,0,nfrom)
  {
    jfrom[i]=readint(fpin);
    jmax=max(jfrom[i],jmax);
    fac[i]=readdouble(fpin);
    printout("normal"," +%lg *coef_c.%d",fac[i],jfrom[i]);
  }
  printout("normal","\n");
  
  noindfwpts=(int *)need("noindfwpts");/* get needed arrays */
  cn=(int *)need("coef_n");
  cc=(double **)need("coef_c");
  nocoefs=(int *)need("nocoefs");
  if (jmax>nocoefs[0]-1)
  { 
    printout("error c_coefadd"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
    exitm4d(0);
  }
  itrange=(int *)find("itrange");
  if (itrange==0) {neqs=0; neqe=noindfwpts[0]; }
  else {neqs=noindfwpts[itrange[0]+1]; neqe=noindfwpts[itrange[1]+2]; }
  
  Loop(ieq,neqs,neqe)
  if (cn[ieq]>0)
  { 
    Loop(i,0,cn[ieq])
    {
      sum=0;
      Loop(j,0,nfrom) sum +=fac[j]*cc[ieq][i+jfrom[j]*cn[ieq]];
      cc[ieq][i+jto*cn[ieq]]=sum;
    }
  }
}

