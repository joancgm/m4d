/* contains c_alias, c_arraylist, arraysize, arraytype, createarray, arraydelete, find, findarrray, need, smalloc, smalloca, tmalloca */

#include "global.h"
#if 0
/* commented out here, real definition in arrays.h */
typedef struct array_
{ 
  char *name;
  int size;        /* size in d,f,i,c or p terms */
  char type;         /* f-float d-double i-int c-char p-pointers   */
  void *pointer;   /*address of array values */
  struct array_ *next;    /* next array =0 if none */
  struct array_ *previous;  /* previous array =0 if none */
  int alias;    /* number of aliases  set to -1 for aliases */
  struct array_ **aliasarray; /* corresponding array pointers */
} Array;
#endif

static Array *array1=0; 
static Array *array_last;
static TArray *tarray_last=0;

/*------------------------------------------------------*/
void c_alias(FILE *fpin, FILE *fprint)
{ 
  char *name; Array *a, *add; int i;
  int *k; double *d; void *p; char *c; float *f; /* consistent with list in smalloca */
  name=readname(fpin);
  printout("normal","alias for %s",name);
  a=findarray(name);
  if (a==0) { printout("warning c_alias",", %s not found, no alias set\n",name); return;}
  if (a->alias==-1) { printout("warning c_alias"," not set %s is an alias\n",name); return;}
  if (a->alias>0) {printout("warning c_alias",", %s already set, no change\n",name); return;}
  a->alias=readint(fpin);
  printout("normal"," %d parts",a->alias);
  if (a->alias <=0) { a->alias=0; return;}
  if ((a->size)%(a->alias)!=0)
  {
    printout("warning c_alias",",%s size %d not divisible by %d, no alias set\n",
             name,a->size,a->alias);
    a->alias=0;
    return;
  }
  a->aliasarray=(Array**)smalloca(a->alias,'p');
  for (i=0;i<a->alias;i++)
  {  
    name=readname(fpin);
    printout("normal"," %s",name);
    add=(Array *)smalloc(sizeof(Array));
    add->name = (char *)smalloc(strlen(name)+1);
    strcpy(add->name,name); 
    a->aliasarray[i]=add;
    add->type=a->type;
    add->size=(a->size)/(a->alias);
    if (a->type=='c') {c=(char *)a->pointer+i*(add->size); add->pointer=(void *)c; }
    if (a->type=='i') {k=(int *)a->pointer+i*(add->size); add->pointer=(void *)k; }
    if (a->type=='d') {d=(double *)a->pointer+i*(add->size); add->pointer=(void *)d; }
    if (a->type=='f') {f=(float *)a->pointer+i*(add->size); add->pointer=(void *)f; }
    if (a->type=='p' || a->type=='s') {p=(void **)a->pointer+i*(add->size); add->pointer=(void **)p; }
    add->alias=-1;
  }
  printout("normal","\n");
}
/*------------------------------------------------------*/
void c_arraylist(FILE *fpin, FILE *fprint)
{ 
  Array *a, *alias; int i;
  a=array1;
  printout("print"," Arrays: name size type pointer\n");
  while (a!=0)
  { 
    printout("print","%-10s %8d  %c  %lu",a->name,a->size,a->type,(long unsigned int)a->pointer);
    if (a->alias>0)
    { 
      printout("print"," alias:");
      for (i=0;i<a->alias;i++)
      {  
        alias=(Array*)(a->aliasarray[i]);
        printout("print"," %s",alias->name);
      }
    }
    printout("print","\n");
    a=a->next;
  }
}
/*---------------------------------------------------------*/
int arraysize(const char *name)
{ 
  int size; Array *array;
  array=findarray(name);
  if (array==0) size=0;
  else size=array->size;
  return (size);
}
/*------------------------------------------------------*/
char arraytype(const char *name)
{ 
  char t; Array *array;
  array=findarray(name);
  if (array==0) t='X';
  else t=array->type;
  return (t);
}
/*------------------------------------------------------*/
void* createarray(const char *name, int size, char type, int nynew)
{  
  Array *newa,*add; void **p,*pointer=0; int i;
  
  newa=findarray(name);  /* check to see if it exists */
  if (newa !=0)
  {
    if (nynew==1 && (size != newa->size || type !=newa->type))
    {
      printout("error createarray","array %s already exists, cannot add it\n",name);     
      exitm4d(0);
    }
    else if (newa->alias==-1)
    { printout("error createarray","%s is being used as an alias, cannot create it\n",name);
      exitm4d(0);
    }
    else
    { 
      if (newa->type=='p' || newa->type=='s') /* free old items in pointer arrays */
      { 
        p=(void **)newa->pointer;
        for (i=0;i<newa->size;i++) if (p[i]!=0) free(p[i]);
      }
      if (size != newa->size || type !=newa->type)
      /* exists with different properties */
      { 
        printout("warning createarray","warning, array %s  replaced",name);
        if (size != newa->size) printout("warning",", size: old %d new %d",
                                       newa->size, size);
        if (type !=newa->type) printout("warning",", type: old %c new %c",
                                      newa->type, type);
        printout("warning","\n");
        free(newa->pointer);
        pointer=smalloca(size,type);
        newa->type=type;
        newa->size=size; 
        newa->pointer=pointer;
      }
      else pointer=newa->pointer;
    }
  }
  
  else /* newa=0  add to list */
  { 
    add=(Array *)smalloc(sizeof(Array));
    add->name = (char *)smalloc(strlen(name)+1);
    strcpy(add->name,name); 
    pointer=smalloca(size,type);
    add->type=type;
    add->size=size; 
    add->pointer=pointer;
    add->alias=0;
    add->next=0;
    add->previous=0;
    if (array1==0) { array1=add; array_last=add;}
    else 
    { 
      array_last->next=(Array *)add;
      add->previous=(Array *)array_last;
      array_last=add;
    } 
  }
  return pointer;
}

int arraydelete(const char *name)
{ 
  Array *a, *ap, *an; int i; char **p;
  a= findarray(name);
  if (a==0) return (0);
  if (a->alias==-1) 
  { 
    printout("warning"," %s is an alias, can only delete primary arrays\n",name);
    return(0);
  }
  if (a->type =='p' || a->type =='s') 
  { 
    p=(char **)a->pointer;
    for (i=0;i<a->size;i++) if (p[i]!=0) free(p[i]);
  }
  free(a->pointer);
  /* check for alias */
  if (a->alias>0) 
  {
    for (i=0;i<a->alias; i++)
    { an=a->aliasarray[i];
      free(an->name);
      free(an);
    }
    free(a->aliasarray);
  }
  ap=a->previous;
  an=a->next;
  if (ap !=0) ap->next=a->next; else array1=an;
  if (an!=0) an->previous=a->previous;  else array_last=ap;
  free(a);
  return(1);
}

/*---------------------------------------------------------*/
void *find(const char *name)
{ 
  void *a; Array *array;
  array=findarray(name);
  if (array==0) a=0;
  else a=array->pointer;
  return (a);
}
/*---------------------------------------------------------*/
Array *findarray(const char *name)
{ 
  Array *find,*alias; int ieq, i;
  if (array1 != 0) 
  {
    find=array1;
    while(find != 0)
    { 
      ieq=strcmp(name,find->name);
      if (ieq == 0) return (find);
      /* check for aliases */
      if (find->alias >0)
        for (i=0;i<find->alias;i++)
        { 
          alias=find->aliasarray[i];
          if (strcmp(name,alias->name) ==0) return (alias);
        }
      find=(Array *)find->next;
    }
  }
  return 0;
}
/*---------------------------------------------------------*/
void *need(const char *name)
{ 
  void *a; Array *array;
  int arrayhowto(const char *name);
  array=findarray(name);
  if (array==0) 
  {
    arrayhowto(name);  /* try arrayhowto */
    
    array=findarray(name);
    if (array==0) 
    {
      printout("error need ","needed array %s not found\n",name); 
      exitm4d(0); 
    }
  }
  a=array->pointer;
  return (a);
}

/*---------------------------------------------------------*/
void *smalloc(int i)    /* safe malloc */
{ 
  void *a;
  if (i<=0)
  { 
    printout("error smalloc","memory request error: %d\n",i);
    exitm4d(0);
  }
  a=malloc(i);
  if (a==NULL)
  { 
    printout("error smalloc","memory allocation error, added memory request: %d\n",i); 
    exitm4d(0);
  }
  return a;
}
/*---------------------------------------------------------*/
void *smalloca(int i, char c)    /* safe malloc  specified type */

{ 
  void *a;
  a=NULL;
  if (c=='c') a=smalloc(i*sizeof(char));
  else if (c=='d') a=smalloc(i*sizeof(double));
  else if (c=='f') a=smalloc(i*sizeof(float));
  else if (c=='i') a=smalloc(i*sizeof(int));
  else if (c=='p' || c=='s') a=smalloc(i*sizeof(void *));
  else 
  {
    printout("error","type %c not listed in smalloca, cannot allocate space\n",c);
    exitm4d(0);
  }
  return a;
}
/*---------------------------------------------------------*/
void *tmalloca(int i, char c)   /* temp malloc, cleared with negative call */

{ 
  void *a; 
  TArray *newa, *before;
  if (i<0) /* clear temp arrays last reacted first deleted */
  {  
    while (tarray_last !=0)
    {
      before=(TArray *)tarray_last->previous;
      free(tarray_last->pointer);
      free(tarray_last);
      tarray_last=before;
      
    }
    return 0;
  }
  a=smalloca(i,c);
  newa=(TArray *)smalloc(sizeof(TArray));
  newa->previous=(void *)tarray_last;
  newa->pointer=a;
  tarray_last=newa;
  return a;
}


