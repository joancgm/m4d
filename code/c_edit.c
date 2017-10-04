/* contains c_edit */
#include "global.h"

void c_edit(FILE *fpin, FILE *fprint)
{
  char *name, type,nameif[10],action, *actionname;
  int nf,i,n,k, idf[10], is[10], ie[10];
  int  iftot, iif[10],ifacf[10], iend, jf;
  Array *af; 
  int *i4d,*mf=0,mv=0; double *df=0,dv=0; char *cf=0,cv=' ',cvold=' ';
  
  
  i4d=(int *)find("idim4d");
  name=readname(fpin);     /* array name */
  type=read1charname(fpin);  /* array type    i d */
  
  nf=readint(fpin);   /* number of dimensions */
  printout("normal"," name %s,   type %c,  %d dims ",name,type,nf);
  if (type !='d' && type !='c' && type !='i')
  { 
    printout("error c_edit","error can't edit array of type %c\n",type); 
    exitm4d(0); 
  }
  iftot=1;
  for (n=0;n<nf;n++)    /* dimension set up */
  { 
    nameif[n]=read1charname(fpin);
    idf[n]=readint(fpin);
    if (i4d>0)
    {
      if (nameif[n]=='i' ) idf[n] += i4d[0];
      else if (nameif[n]=='j') idf[n] += i4d[1];
      else if (nameif[n]=='k') idf[n] += i4d[2];
      else if (nameif[n]=='t') idf[n] += i4d[3];
    }
    iftot*=idf[n];
    if (n==0) ifacf[n]=1;
    else ifacf[n]=ifacf[n-1]*idf[n-1];
    printout("normal"," %c %d",nameif[n],idf[n]);
  }
  printout("normal","\n");
  
  af=findarray(name);
  if (af==0)  /* create array */
  { 
    if (type=='d') 
    {
      df=(double *)createarray(name,iftot,'d',0);
		for (i=0;i<iftot;i++) df[i]=0;
    }
    else if (type=='i') 
    {
      mf=(int *)createarray(name,iftot,'i',0);
		for (i=0;i<iftot;i++) mf[i]=0;
    }
    else if (type=='c') 
    { 
      cf=(char *)createarray(name,iftot,'c',0);
		for (i=0;i<iftot;i++) cf[i]=(char)0;
    }
    printout("normal"," Array %s created and cleared\n",name);
  }
  else
  { 
    if (af->type != type) 
    {
      printout("error c_edit","error, inconsistent type %c %c, can't edit\n",af->type,type); 
      exitm4d(0); 
    }
    if (iftot>af->size) 
    {
      printout("error c_edit","error, size %d larger than array size %d, no edits\n",iftot,af->size); 
      exitm4d(0);
    }
    if (iftot !=af->size) printout("warning"," warning size difference, %d %d\n",iftot,af->size);
    if (af->type=='c') cf=(char *)af->pointer; 
    else if (af->type=='d') df=(double *)af->pointer;
    else if (af->type=='i') mf=(int *)af->pointer;
  }
  /* ok have array, now specify what to do, repeat until action is not valid, 
   all:   s-set;     type i,d:   a-add, m-multiply   type c: reset */
  while(1)
  {  
    actionname=readname(fpin);   
    if (strcmp(actionname,"c:")==0) 
    { 
      fseek(fpin,(long)(-3),1);
      printout("normal","edit: missing end supplied, end edit\n");  
      return; 
    }

    action=actionname[0];
    printout("normal"," %c",action); 
    if (type=='c') 
    { 
      if (action !='s' && action !='r') {printout("normal","  end edit\n"); return;}
      if (action=='r') {cvold=read1charname(fpin); printout("normal","  %c to",cvold); }
      cv=read1charname(fpin); printout("normal","  %c   ",cv);
    }
    else 
    {
      if (action!='s' && action !='a' && action !='m') {printout("normal","  end edit\n"); return;}
      if (type=='d') {dv=readdouble(fpin); printout("normal"," %lg    ",dv);}
      else if (type=='i') {mv=readint(fpin); printout("normal"," %d   ",mv); }
    }
    for (n=0;n<nf;n++)
    {
      is[n]=readint(fpin)-1; if (is[n]<0) is[n]=0; if (is[n]>idf[n]-1) is[n]=idf[n]-1;
      ie[n]=readint(fpin); if (ie[n]>idf[n]) ie[n]=idf[n];
      printout("normal","      %d %d",is[n]+1,ie[n]);
    }
    printout("normal","\n");
    /* set up loop by hand since no of dimensions may vary */
    for (k=0;k<nf;k++) iif[k]=is[k];
    iend=1;
    while(iend)
    { 
      jf=0; for (k=0;k<nf;k++) jf += iif[k]*ifacf[k];
      if (type=='c') 
      { 
        if (action=='s') cf[jf]=cv;
        else if (action=='r' && cf[jf]==cvold) cf[jf]=cv;
      }
      else if (type=='d')
      { 
        if (action=='s') df[jf]=dv;
        else if (action=='a') df[jf]+=dv;
        else if (action=='m') df[jf]*=dv;
      }
      else if (type=='i')
      { 
        if (action=='s') mf[jf]=mv;
        else if (action=='a') mf[jf]+=mv;
        else if (action=='m') mf[jf]*=mv;
      }
      for (k=0;k<nf;k++)
      { 
        iif[k]++;
        if (iif[k]>=ie[k]) iif[k]=is[k]; else break;
        if (k==nf-1 && iif[k]==is[k]) { iend=0; break;}
      }
    }
  }
}

