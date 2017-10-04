/* contains pc_keyword */
/*  analyize an arraydump file to create a lineplot file 
 used for creating convergence plots */

#include "global.h"

void pc_keyword(FILE *fpin, FILE *fprint)
{ 
  char *filein, *fileout;   /* input */
  int nokey; char **keyin; int *nokeyout; char ***keyout;
  
  int  i, j,ktot, n,nline=8, iv;
  char type,cf=' ';
  double a;
  FILE *fin,*fout;
  
  filein=readfilename(fpin);  /* read */
  fileout=readfilename(fpin);  /* read */
  printout("normal","keyword from %s to %s\n",filein,fileout);
  fin=safefopen(filein,"r");
  fout=safefopen(fileout,"w");
  
  nokey=readint(fpin);  /* read */ 
  keyin=(char **)tmalloca(nokey,'p');
  nokeyout=(int *)tmalloca(nokey,'i');
  keyout=(char ***)tmalloca(nokey,'p');
  ktot=0;
  for (i=0;i<nokey;i++)
  {  
    keyin[i]=readname(fpin);  /* read */
    nokeyout[i]=readint(fpin);  /* read */
    ktot+=nokeyout[i];
    keyout[i]=(char **)tmalloca(nokeyout[i],'p');
    for (j=0;j<nokeyout[i];j++) keyout[i][j]=readname(fpin);  /* read */
  }
  
  /* header */
  fprintf(fout," keyword from %s\n %d\n",filein,ktot);
  n=0;
  for (i=0;i<nokey;i++) 
  {
    printout("normal","\n%s:",keyin[i]);
    for (j=0;j<nokeyout[i];j++)
    {
      printout("normal"," %s",keyout[i][j]);
      fprintf(fout," %s",keyout[i][j]);
      n++;
      if (n%nline == 0) printout("normal","\n");
    }
  }
  fprintf(fout,"\n");
  printout("normal","\n");
  ktot=0;
  while (1)
  { 
    n=0;
    for (i=0;i<nokey;i++)
    {
      cf=findword(fin,keyin[i]);
      if (cf==EOF) break;   /* return on end of file */
      iv=readint(fin);
      type=read1charname(fin);
      for (j=0;j<nokeyout[i];j++)
      { 
        a=0;
        if (j<iv) a=readdouble(fin);
        fprintf(fout," %lg",a);
        n++;
        if (n%nline == 0) fprintf(fout,"\n");
      }
    } /* i loop */
    if (cf==EOF) break;
    fprintf(fout,"\n\n");
    ktot++;
  }
  printout("normal"," total %d lines to file %s\n",ktot,fileout);
  fclose(fin);
  fclose(fout);
}

