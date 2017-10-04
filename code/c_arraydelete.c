/* contains c_arraydelete */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_arraydelete(FILE *fpin, FILE *fprint)
{ 
  char *name; int i,j,n;
  n=readint(fpin);
  Loop(j,0,n) 
  { 
    name=readname(fpin);
	 i=arraydelete(name); 
	 if (i==1) printout("normal","array %s deleted\n",name);
	 else printout("normal","array %s not found\n",name);
  }
}
