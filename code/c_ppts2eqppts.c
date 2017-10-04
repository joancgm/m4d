/* contains c_ppts2eqppts */
#include "global.h"
/* take an array defined for equations and set as a full grid points array */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_ppts2eqppts(FILE *fpin, FILE *fprint)
{ 
  int *noindfwpts,*whoisfw; /*needed arrays */
  double *veq;
  double *vg;   /* set */
  char *nameveq,*namevg; /* input */
  int i;
  
  namevg=readname(fpin);
  nameveq=readname(fpin);
  
  noindfwpts=(int *)need("noindppts");
  whoisfw=(int *)need("whoisp");
  vg=(double *)need(namevg);
  veq=(double*)createarray(nameveq,noindfwpts[0],'d',0);
  fprintf(fprint,"%s to %s\n",namevg,nameveq);
  
  Loop(i,0,noindfwpts[0])  veq[i]=vg[whoisfw[i]];
}


