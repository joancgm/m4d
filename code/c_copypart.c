/* contains c_copypart */
#include "global.h"

void c_copypart(FILE *fpin, FILE *fprint)
{
  char *namefrom, *nameto, nameif[10], nameit[10];
  int nf, nt,n,k, idf[10], isf[10], ief[10], idt[10], ist[10], iet[10];
  int idtot, iftot, iif[10], iit[10], ifacf[10], ifact[10], iend, jf, jt;
  Array *af, *at; 
  int *i4d, *dum,*mf=0,*mt=0; float *ff=0,*ft=0; double *df=0,*dt=0; char *cf=0,*ct=0;
  
  
  i4d=(int *)find("idim4d");
  namefrom=readname(fpin);
  nf=readint(fpin);
  printout("normal","from %s  %d\n",namefrom,nf);
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
    ief[n]=min(ief[n],idf[n]);
    isf[n]=min(isf[n],idf[n]);
    if (n==0) ifacf[n]=1;
    else ifacf[n]=ifacf[n-1]*idf[n-1];
    printout("normal"," %c %d %d %d",nameif[n],idf[n],isf[n],ief[n]);
  }
  printout("normal","\n");
  
  af=findarray(namefrom);
  if (af==0) 
  { 
    printout("error c_copypart","error %s does not exist\n",namefrom); 
    exitm4d(0); 
  }
  if (iftot > af->size) 
  {
    printout("error c_copypart","error size larger than array, %d %d\n",iftot,af->size);
    exitm4d(0); 
  }
  
  nameto=readname(fpin);
  nt=readint(fpin);
  printout("normal","to %s %d\n",nameto,nt);
  idtot=1;
  for (n=0;n<nt;n++)
  { 
    nameit[n]=read1charname(fpin);
    idt[n]=readint(fpin);
    ist[n]=readint(fpin);
    iet[n]=readint(fpin);
    if (i4d>0)
    {
      if (nameit[n]=='i' ) idt[n] += i4d[0];
      else if (nameit[n]=='j') idt[n] += i4d[1];
      else if (nameit[n]=='k') idt[n] += i4d[2];
      else if (nameit[n]=='t') idt[n] += i4d[3];
    }
    iet[n]=min(iet[n],idt[n]);
    ist[n]=min(ist[n],idt[n]);
    idtot *= idt[n];
    if (n==0) ifact[n]=1;
    else ifact[n]=ifact[n-1]*idt[n-1];
    printout("normal"," %c %d %d %d",nameit[n],idt[n],ist[n],iet[n]);
  }
  printout("normal","\n");
  at=findarray(nameto);
  if (at==0)
  {      
    printout("normal"," %s not found so it is  being created %c %d\n",nameto,af->type,idtot);
    dum=(int *)createarray(nameto,idtot,af->type,0);
  }
  at=findarray(nameto);
  if (idtot > at->size) 
  {
    printout("error c_copypart","size larger than array, %d %d\n",idtot,at->size);
    exitm4d(0);
  }
  if (at->type != af->type) 
  {
    printout("error c_copypart","error incompatible types %c %d\n",at->type,af->type); 
    exitm4d(0);
  }
  
  if (af->type=='c') {cf=(char *)af->pointer; ct=(char*)at->pointer;}
  else if (af->type=='d') {df=(double *)af->pointer; dt=(double*)at->pointer;}
  else if (af->type=='f') {ff=(float *)af->pointer; ft=(float*)at->pointer;}
  else if (af->type=='i') {mf=(int *)af->pointer; mt=(int *)at->pointer;}
  
  for (k=0;k<nf;k++) iif[k]=isf[k]-1;
  iend=1;
  for (k=0;k<nt;k++) iit[k]=ist[k]-1;
  
  while(iend)
  { 
    jf=0; for (k=0;k<nf;k++) jf += iif[k]*ifacf[k];
    jt=0; for (k=0;k<nt;k++) jt += iit[k]*ifact[k];
    if (af->type=='c') ct[jt]=cf[jf];
    else if (af->type=='d') dt[jt]=df[jf];
    else if (af->type=='f') ft[jt]=ff[jf];
    else if (af->type=='i') mt[jt]=mf[jf];
    for (k=0;k<nt;k++)
    { 
      iit[k]++;
      for (n=0;n<nf;n++)
      {
        if (nameif[n]==nameit[k])
        { iif[n]++;
          if (iif[n]==ief[n]) iif[n]=isf[n]-1;
        }
      }
      if (iit[k]==iet[k]) iit[k]=ist[k]-1; else break;
      if (k==nt-1 && iit[k]==ist[k]-1) iend=0;
    }
    if (iend==0) break;
  }
}
