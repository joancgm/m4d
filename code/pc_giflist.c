/* contains pc_giflist */
#include "global.h"
#include "p_plot.h"
/* create a list of gif files prename istart iend */

/* ---------------------------- */
void pc_giflist(FILE *fpin, FILE *fprint)
{
  int i,is,ie;
  char *prename,*listname;
  FILE *fp;
  
  prename=readname(fpin);
  is=readint(fpin);
  ie=readint(fpin);
  printout("normal"," %s %d %d\n",prename,is,ie);
  listname=setname2(prename,".list");
  fp=fopen(listname,"w");
  for(i=is;i<=ie;i++) fprintf(fp,"%s%d.gif\n",prename,i);
  fclose(fp);
}
