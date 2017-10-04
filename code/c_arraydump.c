/* contains c_arraydump, c_arraydumpmore arraydump */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

void arraydump(FILE *fpin, FILE *fprint, FILE *fp);

/* dump arrays to named file */
void c_arraydump(FILE *fpin, FILE *fprint)
{ 
  char *fileout;
  FILE *fp;
  
  fileout=readfilename(fpin);
  fp=safefopen(fileout,"w");
  printout("normal","dumping to file %s:  ",fileout);
  arraydump(fpin,fprint,fp);
  fclose(fp);	
}
/* dump append arrays to named file */
void c_arraydumpmore(FILE *fpin, FILE *fprint)
{ 
  char *fileout;
  FILE *fp;
  
  fileout=readfilename(fpin);
  fp=safefopen(fileout,"a");
  printout("normal","appending dump to file %s:  ",fileout);
  arraydump(fpin,fprint,fp);
  fclose(fp);
}
void arraydump(FILE *fpin, FILE *fprint, FILE *fp)
{  
  int i,n,*iv;
  char *name,*cv;
  double *v;
  Array *array;
  
  n=100;
  while (n)
  { 
    name=readname(fpin);
    if (name[0]=='\0')  {printout("normal","\n");  break; }
    if (strcmp(name,"c:")==0) { fseek(fpin,(long)(-3),1);printout("normal","\n");  break; }
    n--;
    array=findarray(name);
    if (array==0)  {printout("warning arraydump"," array %s not found\n",name); continue; }
    printout("normal"," %s",name);
    if (array->type=='p') {printout("warning arraydump"," can't (yet) dump pointer arrays\n"); continue; }
    fprintf(fp,"%s %d %c",array->name, array->size, array->type);
    if (array->type=='d')
    { 
      v=(double *)array->pointer;
      Loop(i,0,array->size)
      { 
        if ((i%10)==0 && array->size > 1) fprintf(fp,"\n");
        fprintf(fp," %lg",v[i]); 
      }
      fprintf(fp,"\n");
    }
    else if (array->type=='i')
    { 
      iv=(int *)array->pointer;
      Loop(i,0,array->size)
      { 
        if ((i%20)==0 && array->size > 1) fprintf(fp,"\n");
        fprintf(fp," %d",iv[i]); 
      }
      fprintf(fp,"\n"); 
    }
    else if (array->type=='c')
    { 
      cv=(char *)array->pointer;
      Loop(i,0,array->size)
      { 
        if ((i%20)==0 && array->size > 1) fprintf(fp,"\n");
        fprintf(fp," %c",cv[i]); 
      }
      fprintf(fp,"\n"); 
    }
  }
}
