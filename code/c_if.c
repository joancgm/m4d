/* contains c_if */
#include "global.h"

void c_if(FILE *fpin, FILE *fprint)
{ 
  char *nif[2],*leg,*ifpass,*stringp=0,*iffail,*stringf=0;  /* read */ 
  
  double diff[2],*dif[2]; int i,*iff[2],nypass,iiff[2],j;
  char type[2],c;
  /* if nif1[0] <=>%(leg) nif2[0], 
   is true either find a string(ifpass="find") or continue, 
   if fail either find a string (iffail="find") or continue
   if nif[0] or nif[1] does not exist or leg is not valid, the fail option is used
   */
  
  nif[0]=readname(fpin);
  leg=readname(fpin);
  nif[1]=readname(fpin);
  ifpass=readname(fpin);
  printout("normal"," if %s %s %s, then %s",nif[0],leg,nif[1],ifpass);
  if (strcmp(ifpass,"find")==0) 
  { stringp=readname(fpin); printout("normal"," %s",stringp); }
  iffail=readname(fpin);
  printout("normal",", otherwise %s",iffail);
  if (strcmp(iffail,"find")==0) 
  { stringf=readname(fpin); printout("normal"," %s",stringf); }
  
  
  type[0]=arraytype(nif[0]);
  type[1]=arraytype(nif[1]);
  nypass=0;
  for (i=0;i<2;i++)
  { 
    if (type[i]=='d') 
    { dif[i]=(double *)need(nif[i]); diff[i]=dif[i][0]; }
    else if (type[i]=='i')
    { iff[i]=(int *)need(nif[i]); iiff[i]=iff[i][0]; diff[i]=iiff[i]; }
    else if (type[i]=='X') /* try to interpret at value */
    { 
      j=sscanf(nif[i],"%d",&iiff[i]);
      j+=sscanf(nif[i],"%lg",&diff[i]);
      if (j==0)
        {
          printout("normal","  Failed: %s not found\n",nif[i]); 
          nypass=-1;
        }
    }
  }
  if (nypass==0)
  {
    if (type[0]=='i' || type[1]=='i')
    { 
      printout("normal"," %d %c %d ? ",iiff[0],leg[0],iiff[1]);
      if (leg[0]=='<') {if (iiff[0]<iiff[1]) nypass=1; }
      else if (leg[0]=='=') {if (iiff[0]==iiff[1]) nypass=1;}
      else if (leg[0]=='>') {if (iiff[0]>iiff[1]) nypass=1; }
      else if (leg[0]=='%') {if (iiff[0]%iiff[1]==0) nypass=1;}
      else {printout("warning c_if","  Failed: operation %s not valid\n",leg); nypass=-2;} 
    }
    else 
    { 
      printout("normal","% lg %c %lg ? ",diff[0],leg[0],diff[1]);
      if (leg[0]=='<') {if (diff[0]<diff[1]) nypass=1; }
      else if (leg[0]=='=') {if (diff[0]==diff[1]) nypass=1;}
      else if (leg[0]=='>') {if (diff[0]>diff[1]) nypass=1; }
      else {printout("warning","  Failed: operation %s not valid for type 'd'\n",leg); nypass=-2;} 
    }
  }
  /*if (nypass==-2) { printout("normal"," exit on bad input\n"); exitm4d(0); } */
  
  if (nypass==1)
  {
    if (strcmp(ifpass,"find")==0) 
    { 
      c=findstring(fpin,stringp);
      if (c==EOF) printout("warning c_if"," True but warning %s not found, EOF\n",stringp);
      else printout("normal"," True, %s found\n",stringp);
    }
    else printout("normal"," True, %s\n",ifpass);
    return;
  }
  else 
  { 
    if (strcmp(iffail,"find")==0) 
    { 
      c=findstring(fpin,stringf);
      if (c==EOF) printout("warning c_if"," Not true but warning %s not found, EOF\n",stringf);
      else printout("normal"," Not true, %s found\n",stringf);
    }
    else printout("normal"," Not true, %s\n",iffail);
    return;
  }
}
