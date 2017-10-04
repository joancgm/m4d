/* contains c_gridtomefp */
#include "global.h"

/* dump  mefp style grid for plotting  read file name  then  read t index  */
void c_gridtomefp(FILE *fpin, FILE *fprint)
{   
  double *xyz, *abcd,*aa;
  int *i4d, it,i,j,k,iall,m;
  char *fileout,*clt,cc;
  FILE *fp;
  
  fileout=readfilename(fpin);
  fp=safefopen(fileout,"w");
  it=readint(fpin);
  xyz=(double *)need("xyz");
  abcd=(double *)need("abcd");
  i4d=(int *)need("idim4d");
  clt=(char *)need("clt");
  it=min(it,i4d[3]);
  iall=i4d[0]*i4d[1]*i4d[2];
  fprintf(fp,"grid at it=%d\n%d %d %d %d\n",it,i4d[0],i4d[1],i4d[2],1);
  aa=abcd+i4d[0];  /* B */
  for (i=0;i<i4d[1];i++)
  {fprintf(fp," %lg ",(float)aa[i]); if ((i+1)%10==0) fprintf(fp,"\n");}
  fprintf(fp,"\n");	
  aa+=i4d[1];  /* C */
  for (i=0;i<i4d[2];i++)
  {fprintf(fp," %lg ",(float)aa[i]); if ((i+1)%10==0) fprintf(fp,"\n");}
  fprintf(fp,"\n");	
  aa=abcd;
  for (i=0;i<i4d[0];i++)  /* each a plane */
  {
    fprintf(fp," %lg \n",(float)aa[i]); 
    for (k=0;k<i4d[2];k++) for (j=0;j<i4d[1];j++)
    {
      m=i+i4d[0]*(j+i4d[1]*(k+i4d[2]*it));
      cc=clt[m];	
      if (cc=='f') fprintf(fp," 1");
      else if (cc=='w') fprintf(fp," 3");
      else if (cc=='s') fprintf(fp," 5");
      else fprintf(fp," 2");
      fprintf(fp," %lg %lg %lg\n",xyz[m],xyz[m+iall*i4d[3]],xyz[m+2*iall*i4d[3]]);
    }
  }
  printout("normal"," 3d grid at it= %d dumped mefp style to file %s\n",it,fileout);
  fclose(fp);
}
