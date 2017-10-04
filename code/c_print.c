/* contains c_print */
#include "global.h"

/* print arrays */
void c_print(FILE *fpin, FILE *fprint)
{   
  char *name, cdfi,*cblock=0, **sblock=0;
  int *iblock=0,norep,isize,*i4d,iall,iformat;
  double *dblock=0; float *fblock=0;
  int i,i1,i2,i3,irep,jd[4],m;
  Array *array;
  
  name=readname(fpin);
  iformat=readint(fpin);
  
  array=findarray(name);
  if (array==0)  {printout("print","print: array %s not found\n",name); return;}
  
  cdfi=array->type;
  if (cdfi=='f') fblock=(float *)array->pointer;
  else if (cdfi=='d') dblock=(double *)array->pointer;
  else if (cdfi=='c') cblock=(char*)array->pointer;
  else if (cdfi=='i') iblock=(int *)array->pointer;
  else if (cdfi=='s') sblock=(char **)array->pointer;
  else {printout("print","print: type %c not set up for print\n",cdfi); return;}
  
  isize=array->size;
  printout("print","\n  %s  size %d  type %c",name,isize,cdfi);
  if (isize<11 || cdfi=='s')  /* small array or string */
  {
    printout("print"," values:");
    if (cdfi=='d')  for (m=0;m<isize;m++) printout("print"," %lg",dblock[m]);
    else if (cdfi=='i') for (m=0;m<isize;m++) printout("print"," %d",iblock[m]);
    else if (cdfi=='c') for (m=0;m<isize;m++) printout("print"," %c",cblock[m]);
    else if (cdfi=='s') for (m=0;m<isize;m++) printout("print"," %s",sblock[m]);
    printout("print","\n");
    return;
  }
  printout("print","\n");
  i4d=(int *)need("idim4d");
  for (i=0;i<4;i++) jd[i]=i4d[i]+iformat;
  if (iformat==2) for (i=0;i<3;i++) jd[i]=2*i4d[i]-1;
  jd[3]=i4d[3];
  iall=jd[0]*jd[1]*jd[2]*jd[3];
  if (iformat<=2 && isize%iall!=0) iformat=10;
  
  if (iformat>2)   
  {
    printout("print","array not standard or specified size\n"); 
    for (m=0;m<isize;m++) 
    {
      if (cdfi=='f') printout("print"," %g",fblock[m]);
      else if (cdfi=='d') printout("print"," %lg",dblock[m]);
      else if (cdfi=='i') printout("print"," %d",iblock[m]);
      else if (cdfi=='c') printout("print"," %c",cblock[m]);
      if ((m+1)%abs(iformat)==0) printout("print","\n");
    }
    printout("print","\n");
    return;
  }
  
  norep=isize/iall;
  for (irep=0;irep<norep;irep++)
  {
    printout("print","\n%s repeat  %d dimensions %d %d %d %d\n",
            name,irep,jd[0],jd[1],jd[2],jd[3]);
    for (i3=0;i3<jd[3];i3++) for (i2=0;i2<jd[2];i2++) for (i1=0;i1<jd[1];i1++)
    {
      printout("print","%4d %4d %4d    ",i1,i2,i3);
      for (i=0;i<jd[0];i++)
      {
        m=irep*iall+i+jd[0]*(i1+jd[1]*(i2+jd[2]*i3));
        if (cdfi=='f') printout("print"," %g",fblock[m]);
        else if (cdfi=='d') printout("print"," %lg",dblock[m]);
        else if (cdfi=='i') printout("print"," %d",iblock[m]);
        else if (cdfi=='c') printout("print"," %c",cblock[m]);
      }
      printout("print","\n");
    }
  }
}
