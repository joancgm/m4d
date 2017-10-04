/* contains c_prints */
#include "global.h"

void c_prints(FILE *fpin, FILE *fprint)
{
  char *name, nameif[10];
  int nf,n,k, idf[10], isf[10], ief[10];
  int  iftot, iif[10],ifacf[10], iend, jf;
  Array *af; 
  int *i4d,*mf=0; float *ff=0; double *df=0; char *cf=0;
  
  
  i4d=(int *)find("idim4d");
  name=readname(fpin);
  nf=readint(fpin);
  printout("print"," %s  %d ",name,nf);
  iftot=1;
  for (n=0;n<nf;n++)
  { 
    nameif[n]=read1charname(fpin);
    idf[n]=readint(fpin);
    isf[n]=readint(fpin);
    ief[n]=readint(fpin);
    if (i4d>0)
    {
      if (nameif[n]=='i' ) idf[n] += i4d[0];
      else if (nameif[n]=='j') idf[n] += i4d[1];
      else if (nameif[n]=='k') idf[n] += i4d[2];
      else if (nameif[n]=='t') idf[n] += i4d[3];
    }
    iftot*=idf[n];
    ief[n]=max(ief[n],isf[n]);
    ief[n]=min(ief[n],idf[n]);
    isf[n]=min(isf[n],idf[n]);
    if (n==0) ifacf[n]=1;
    else ifacf[n]=ifacf[n-1]*idf[n-1];
    printout("print"," %c %d %d %d",nameif[n],idf[n],isf[n],ief[n]);
  }
  
  af=findarray(name);
  if (af==0) {printout("print","   %s not found\n",name); return; }
  printout("print"," type  %c\n",  af->type);
  if (iftot !=af->size) printout("print"," warning size difference, %d %d\n",iftot,af->size);
  if (iftot>af->size) {printout("print"," not printing\n"); return; }
  if (af->type=='c') cf=(char *)af->pointer; 
  else if (af->type=='d') df=(double *)af->pointer;
  else if (af->type=='f') ff=(float *)af->pointer; 
  else if (af->type=='i') mf=(int *)af->pointer;
  else { printout("print"," print not implemented for type %c\n",af->type); return;}
  printout("print"," across %c from %d   through %d\n",nameif[0], isf[0]-1,ief[0]-1);
  for (k=0;k<nf;k++) iif[k]=isf[k]-1;
  iend=1;
  
  while(iend)
  {
    jf=0; for (k=0;k<nf;k++) jf += iif[k]*ifacf[k];
    if (iif[0]==isf[0]-1) 
    { 
      printout("print","\n"); 
      for (k=1;k<nf;k++) printout("print"," %d",iif[k]); 
    }
    if (af->type=='c') printout("print"," %c",cf[jf]);
    else if (af->type=='d') printout("print"," %lg",df[jf]);
    else if (af->type=='f') printout("print"," %g",ff[jf]);
    else if (af->type=='i') printout("print"," %d",mf[jf]);
    for (k=0;k<nf;k++)
    {
      iif[k]++;
      if (iif[k]==ief[k]) iif[k]=isf[k]-1; else break;
      if (k==nf-1 && iif[k]==isf[k]-1) { iend=0; break;}
    }
  }
  printout("print","\n");
}
