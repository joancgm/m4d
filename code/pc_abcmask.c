/* contains pc_abcmask */
#include "global.h"
#include "p_plot.h"
/*  set points in an array name to value for specified a, b, or c values
 input: name, value, aborc, abcvalues until value less than previous
 */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* ---------------------------- */
void pc_abcmask(FILE *fpin, FILE *fprint)
{
  int *i4d; double *a[4],*prop; /* need */
  char *name,*aborc;
  double value,abcvalue[500],*abc,f;
  int i,L,L2,L3,ii[4],iv,ivmax=0;
  
  /* get needed arrays */
  i4d=(int *)need("idim4d");
  a[0]=(double *)need("abcd");
  for (i=1;i<4;i++) a[i]=a[i-1]+i4d[i-1];
  
  /* read input */
  name=readname(fpin);
  value=readdouble(fpin);
  aborc=readname(fpin);
  printout("normal","set %s = %lg for %s = \n",name,value,aborc);
  
  for (i=0;i<500;i++)
  {
    abcvalue[i]=readdouble(fpin);
    if (i>0 && abcvalue[i]<abcvalue[i-1]) {ivmax=i; break; }
    printout("normal"," %lg",abcvalue[i]);
  }
  printout("normal","\n");
  prop=(double *)need(name);
  if (strcmp(aborc,"a") ==0) {L=0; abc=a[0];}
  else if (strcmp(aborc,"b") ==0) {L=1; abc=a[1];}
  else if (strcmp(aborc,"c") ==0) {L=2; abc=a[2];}
  else
  {
    abc=(double*)need(aborc);
    if (aborc[0]=='a') L=0;
    else if (aborc[0]=='b') L=1;
    else if (aborc[0]=='c') L=2;
    else
    {
      printout("warning pc_abcmask"," warning:  %s does not begin with a, b, or c, no action\n",aborc);
      return;
    }
  }
  printout("normal","for ii[%d]=\n",L);
  L2=(L+1)%3;
  L3=(L+2)%3;
  ii[3]=0;
  for (iv=0;iv<ivmax;iv++)
  {
    ii[L]=findex(abcvalue[iv],abc,i4d[L],&f);
    if (f>.5) ii[L]++;
    printout("normal"," %d",ii[L]);
    for (ii[L2]=0;ii[L2]<i4d[L2];ii[L2]++)
      for (ii[L3]=0;ii[L3]<i4d[L3];ii[L3]++)
        prop[In4(ii,i4d)]=value;
  }
  printout("normal","\n");
}