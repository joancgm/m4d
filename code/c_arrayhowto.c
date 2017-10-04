/* contains c_arrayhowto arrayhowto */

#include "global.h"
#include "commandlist.h"

static char *arrayhowtofile=0;

/*------------------------------------------------------*/
void c_arrayhowto(FILE *fpin, FILE *fprint)
{  
  char *name; FILE *fpinh;
  name=readfilename(fpin);
  fpinh=fopen(name,"r");
  if (fpinh==NULL) printout("warning c_arrayhowto",
                            "file %s not found, no arrayhowtofile set\n",name);
  else
  { 
    arrayhowtofile=(char *)smalloca(strlen(name)+1,'c');
    strcpy(arrayhowtofile,name);
    printout("normal","arrayhowtofile set to %s\n",arrayhowtofile);
    fclose(fpinh);
  }
}
/* ----------------------- */
int arrayhowto(const char *name)
{ 
  FILE *fpinh;
  char c,*nn,*com;
  int i=0;
  char **info4print;
  
  if (arrayhowtofile==0) return 0;
  fpinh=fopen(arrayhowtofile,"r");
  if (fpinh==NULL) return 0;
  while(1)
  { 
    c=findstring(fpinh,"v:");
    if (c==EOF) return 0;
    nn=readname(fpinh);
    if (strcmp(nn,name)!=0) continue;
    com=readname(fpinh);
    printout("normal","****  arrayhowto making %s using %s\n",name,com);
    info4print=(char **)need("info4print");
    i=commandlist(com,fpinh,(FILE *)info4print[2]);
    if (i==0) printout("normal","        ***** FAILED %s not found|n",com);
    else printout("normal"," OK\n");
    fclose(fpinh);
    return i;
  }
}
