/* contains c_varmatch */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

/* set variables so consistent with matched points */
void c_varmatch(FILE *fpin, FILE *fprint)
{ 
  int *npts,*match; double *ch; char *name;
  int i,j,n,size,irep;
  npts=(int *)need("nogpts");   /* get needed arrays */
  match=(int *)need("match");
  n=100;
  while (n)
  {
    name=readname(fpin);
    n--;
    if (name[0]=='\0') break;
    if (strcmp(name,"c:")==0) { fseek(fpin,(long)(-3),1); break; }
    ch=(double *)find(name);
    if (ch==0) printout("warning c_varmatch","%s not found\n",name);
    else
    { 
      size=arraysize(name);
      if (size%npts[0] !=0) 
      {
        printout("warning c_varmatch"," wrong size %d, numpts %d\n",size,npts[0]);
        continue;
      }
      irep=size/npts[0];
      printout("normal","matching %s\n",name);
      Loop(j,0,irep)
      Loop(i,0,npts[0]) ch[i+j*npts[0]]=ch[match[i]+j*npts[0]];
    }
  }
}
