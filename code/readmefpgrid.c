/* contains readmefpgrid */
#include "global.h"

/* read a 3-d MEFP style grid if lt=2 set to 1 if lt=4 set to 3 */

void readmefpgrid(char *filename, double **axyz, int **alt, double **aabc, int *idim)
{
  char *commentline; 
  double *a, *b, *c, *xyz,*abc; 
  int i,j,k,m,L,iall,id,jd,kd,lcoor,isum,*lt;
  FILE *fp;
  
  fp=safefopen(filename,"r");
  printout("normal"," master geometry: MEFP style gridfile= %s\n",filename);
  
  /*---------------comment line------------*/ 
  commentline=readline(fp);
  printout("normal","%s\n",commentline);
  id=readint(fp);
  jd=readint(fp);
  kd=readint(fp);
  lcoor=readint(fp);
  printout("normal","%d %d %d %d\n",id, jd, kd, lcoor);
  iall=id*jd*kd;
  isum=id+jd+kd;
  xyz=(double *)tmalloca(iall*3,'d');
  lt=(int *)tmalloca(iall,'i');
  abc=(double *)tmalloca(isum,'d');
  a=abc; b=a+id; c=b+jd;
  *axyz=xyz; *alt=lt; *aabc=abc; 
  idim[0]=id; idim[1]=jd; idim[2]=kd;
  
  for (j=0;j<jd;j++) b[j]=readdouble(fp);
  for (k=0;k<kd;k++) c[k]=readdouble(fp);
  
  for (i=0;i<id;i++)
  {
    a[i]=readdouble(fp);
    for (k=0;k<kd;k++) for (j=0;j<jd;j++)
    {
      m=i+id*(j+jd*k);
		lt[m]=readint(fp);
      if (lt[m]==2) lt[m]=1;
      if (lt[m]==4) lt[m]=3;
      for (L=0;L<3;L++) xyz[m+iall*L]=readdouble(fp);
    }
  }
  printout("normal","a= %lg to %lg, b= %lg to %lg, c= %lg to %lg\n",a[0],a[id-1],
         b[0],b[jd-1],c[0],c[kd-1]); 
  fclose(fp);
}

