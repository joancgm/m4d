/* contains c_coefcplus */
#include "global.h"
/* set namecplus as sum of positive coefficients is set jcoef  */

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_coefcplus(FILE *fpin, FILE *fprint)
{ 
  int *noindfwpts,*cn,*nocoefs;/* needed arrays */
  double **cc; 
  int *itrange;  /* use if available */
  double *cplus;  /* set */
  int jcoef; char *namecplus;  /* input parameters */
  
  int i,ieq;
  int neqs,neqe;
  
  jcoef=readint(fpin);
  namecplus=readname(fpin);
  printout("normal","set %s sum of + coef set %d\n",namecplus,jcoef);
  noindfwpts=(int *)need("noindfwpts"); /* get needed arrays */
  nocoefs=(int *)need("nocoefs");
  if (jcoef>nocoefs[0]-1)
  { 
    printout("error c_coefcplus"," error, coef only dimensioned for %d coefs\n",nocoefs[0]);
    exitm4d(0);
  }
  cn=(int *)need("coef_n");
  cc=(double **)need("coef_c");
  itrange=(int *)find("itrange");
  if (itrange==0) {neqs=0; neqe=noindfwpts[0]; }
  else {neqs=noindfwpts[itrange[0]+1]; neqe=noindfwpts[itrange[1]+2]; }
  
  cplus=(double *)createarray(namecplus,noindfwpts[0],'d',0);
  
  Loop(ieq,neqs,neqe)
  { 
    cplus[ieq]=0;
    if (cn[ieq]>0)
		Loop(i,0,cn[ieq])
      if (cc[ieq][i+jcoef*cn[ieq]]>0) 
        cplus[ieq]+=cc[ieq][i+jcoef*cn[ieq]];
  }
}



