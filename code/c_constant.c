/* contains c_constant */
#include "global.h"
/* set constants or small arrays of constants */
#define Loop(n,a,b) for (n=a;n<b;n++)

void c_constant(FILE *fpin, FILE *fprint)
{ 
  char *name,type; int length;
  int i,*iv; double *vv; char *cv,**sv;
  
  name=readname(fpin);
  type=read1charname(fpin);
  length=readint(fpin);
  printout("normal","%s, type: %c, length %d, values",name,type,length);
  if (type=='i')
  {
    iv=(int *)createarray(name,length,'i',0);
	 Loop(i,0,length)   { iv[i]=readint(fpin); printout("normal"," %d",iv[i]); }
  }
  else if (type=='d')
  {
    vv=(double *)createarray(name,length,'d',0);
	 Loop(i,0,length)   { vv[i]=readdouble(fpin); printout("normal"," %lg",vv[i]); }
  }
  else if (type=='c')
  {
    cv=(char *)createarray(name,length,'c',0);
	 Loop(i,0,length)   { cv[i]=read1charname(fpin); printout("normal"," %c",cv[i]); }
  }
  else if (type=='s')
  {
    sv=(char **)createarray(name,length,'s',0);
	 Loop(i,0,length)   
    { 
      sv[i]=readnamep(fpin); 
      printout("normal"," %s",sv[i]); 
    }
  }
  else printout("warning c_constant"," allowed types are i(int) d(double) c(char) s(char *)");
  printout("normal","\n"); 
}
