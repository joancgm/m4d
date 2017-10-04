/* contains c_eqppts2ppts */
#include "global.h"
/* take an array defined for p-equations and set as a full p-points array */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_eqppts2ppts(FILE *fpin, FILE *fprint)
{ 
  int *noindppts,*whoisp, *noppts, *matchpc=0; /*needed arrays */
  double *cpsleep=0;
  double *veq;
  double *vg;   /* set */
  char *nameveq,*namevg; /* input */
  int nymatch;
  char *nyset;
  
  int i,j,k,n=0,isize,irep,m1,m2;
  
  nameveq=readname(fpin);
  namevg=readname(fpin);
  nymatch=readint(fpin);
  
  noindppts=(int *)need("noindppts");
  noppts=(int *)need("noppts");
  whoisp=(int *)need("whoisp");
  matchpc=(int *)need("matchpc");
  if (nymatch>1) cpsleep=(double *)need("cpsleep");
  veq=(double *)need(nameveq);
  isize=arraysize(nameveq);
  irep=1;
  if (isize%noindppts[0]==0) irep=isize/noindppts[0];
  else 
  { 
    printout("warning c_eqppts2ppts","size error for %s ignored size %d not divisible by %d, no matching\n",
            nameveq,isize,noindppts[0]);
    nymatch=0;
  }
  vg=(double*)createarray(namevg,noppts[0]*irep,'d',0);
  printout("normal","%s to %s,  irep=%d, nymatch %d\n",nameveq,namevg,irep,nymatch);
  
  nyset=(char *)tmalloca(noppts[0],'c');
  Loop(i,0,noppts[0]) nyset[i]='n';
  Loop(i,0,noppts[0]*irep) vg[i]=0;
  Loop(i,0,noindppts[0])
  {  
    nyset[whoisp[i]]='y';
    Loop(j,0,irep) 
    vg[whoisp[i]+j*noppts[0]]=veq[i+j*noindppts[0]];
  }
  
  if (nymatch>0) 
    Loop(k,0,10)
  { 
    n=0;
    Loop(i,0,noppts[0]) 
    {
      if (nyset[i]=='y') continue;
      n++;
      m1=matchpc[i];
      if (m1<0 || nyset[m1]=='n') continue;
      m2=matchpc[i+noppts[0]];
      if (nyset[m2]=='n') continue;
      n--;
      nyset[i]='y';
      Loop(j,0,irep) 
      { 
        if (nymatch==1) vg[i+j*noppts[0]]=.5*(vg[m1+j*noppts[0]]+vg[m2+j*noppts[0]]);
        else vg[i+j*noppts[0]]=cpsleep[i]*vg[m1+j*noppts[0]]+cpsleep[i+noppts[0]]*vg[m2+j*noppts[0]];
      }
    }
    if (n==0) break;
  }
  if (n>0) printout("normal"," %d sleeping points unresolved\n",n);
}


