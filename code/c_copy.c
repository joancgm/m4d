/* contains c_copy */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_copy(FILE *fpin, FILE *fprint)
{ 
  char *namefrom, *nameto;  
  Array *af; double*v,*vf; int i,*iv,*ivf; char *cv,*cvf;
  
  namefrom=readname(fpin);
  nameto=readname(fpin);
  printout("normal","copy %s to %s",namefrom,nameto);
  af=findarray(namefrom);
  if (af==0) 
  {
    printout("exit c_copy",": cannot copy, %s does not exist\n",namefrom);
    exitm4d(0);
  }
  else
  { 
    printout("normal","  size=%d, type=%c\n",af->size,af->type);
	 if (af->type=='d') 
    { 
      v=(double *)createarray(nameto,af->size,af->type,0);
      vf=(double *)af->pointer;
      Loop(i,0,af->size) v[i]=vf[i];
    }
	 else if (af->type=='i') 
    { 
      iv=(int *)createarray(nameto,af->size,af->type,0);
      ivf=(int *)af->pointer;
      Loop(i,0,af->size) iv[i]=ivf[i];
    }
	 else if (af->type=='c') 
    { 
      cv=(char *)createarray(nameto,af->size,af->type,0);
      cvf=(char *)af->pointer;
      Loop(i,0,af->size) cv[i]=cvf[i];
    }
	 else 
    {
      printout("error c_copy",": cannot copy array of type %c\n",af->type);
      exitm4d(0);
    }
  }
}
