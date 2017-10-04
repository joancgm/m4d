/* contains c_arrayread */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_arrayread(FILE *fpin, FILE *fprint)
{ 
  int i,*iv,nyneed,size,ok;
  char *filein,*name,*cv,type;
  double *v;
  FILE *fp;
  
  nyneed=readint(fpin);
  filein=readfilename(fpin);
  fp=fopen(filein,"r");
  if (fp==NULL) 
  { 
    printout("warning c_arrayread"," file %s not found",filein);
    if (nyneed>0) 
    {
      printout("error c_arrayead",",   nyneed %d STOP\n",nyneed); 
      exitm4d(0); 
    }
    printout("warning","\n");
    return;
  }
  printout("normal"," reading from file %s\n  ",filein);
  ok=1;
  while(ok)
  { 
    name=readname(fp);
    if (name==NULL) break;
    if (name[0]=='\0') break;
    size=-1;
    type='x';
    size=readint(fp);
    type=read1charname(fp);
    printout("normal","           %s   %d   %c\n",name,size,type);
    if (size<1) {ok=0; break; }
    if (type != 'd' && type !='c' && type !='i') {ok=0; break; }
    if (type=='d') 
    { 
      v=(double *)createarray(name,size,type,0);
      Loop(i,0,size) {v[i]=0; v[i]=readdouble(fp); }
    }
    else if (type=='i') 
    { 
      iv=(int *)createarray(name,size,type,0);
      Loop(i,0,size) { iv[i]=0; iv[i]=readint(fp); }
    }
    else if (type=='c') 
    { 
      cv=(char *)createarray(name,size,type,0);
      Loop(i,0,size) { cv[i]='_'; cv[i]=read1charname(fp); }
    }
  }
  if (ok==0)
  { 
    printout("warning c_arrayread"," input error\n");
    if (nyneed>0) 
    {
      printout("error c_arrayread","STOP\n"); 
      exitm4d(0); 
    }
  }
  fclose(fp);
}
