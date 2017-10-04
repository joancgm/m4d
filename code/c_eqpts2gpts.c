/* contains c_eqpts2gpts */
#include "global.h"
/* take an array defined for equations and set as a full grid points array */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_eqpts2gpts(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*noindfwpts,*whoisfw, *match; /*needed arrays */
  double *veq;
  double *vg;   /* set */
  char *nameveq,*namevg; /* input */
  int nymatch;
  
  int i,j,iall,isize,irep;
  
  nameveq=readname(fpin);
  namevg=readname(fpin);
  nymatch=readint(fpin);
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  noindfwpts=(int *)need("noindfwpts");
  match=(int *)need("match");
  whoisfw=(int *)need("whoisfw");
  veq=(double *)need(nameveq);
  isize=arraysize(nameveq);
  irep=1;
  if (isize%noindfwpts[0]==0) irep=isize/noindfwpts[0];
  else printout("warning c_eqpts2gpts","size error for %s ignored\n",nameveq);
  vg=(double*)createarray(namevg,iall*irep,'d',0);
  printout("normal","%s to %s, nymatch=%d  irep=%d\n",nameveq,namevg,nymatch,irep);
  
  Loop(i,0,iall*irep) vg[i]=0;
  Loop(i,0,noindfwpts[0]) 
  Loop(j,0,irep) 
  vg[whoisfw[i]+j*iall]=veq[i+j*noindfwpts[0]];
  if (nymatch>0) 
    Loop(i,0,iall) 
    Loop(j,0,irep) vg[i+j*iall]=vg[match[i]+j*iall];
  
}


