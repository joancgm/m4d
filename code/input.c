/* contains findstring, findword, readdouble, readfilename, readfloat, readint, readline, readname, readnamep, read1charname, safefopen, setname, setname2, setname3 */
#include "global.h"

/*---------------------------------------------------------*/
char findstring(FILE *fpin, const char *st)
{ 
  char c='x';
  int  i,j,k;
  
  i=strlen(st);
  for (j=0;j<i;k=1)
  {
    c=getc(fpin);
    if (c==EOF) j=i;
    else if (c==st[j]) j++;
    else j=0;
  }
  return c;
}/*---------------------------------------------------------*/
char findword(FILE *fpin, const char *st)
{ 
  char format[20], *string;
  int  i,j=EOF;
  
  i=strlen(st);
  string=smalloca(i+2,'c');
  string[0]='\0';
  sprintf(format,"%%%ds",i+1);
  while (strcmp(st,string)!=0)
  {
    j=fscanf(fpin,format,string);
    if (j==EOF) break;
  }
  free(string);
  if (j==EOF) return EOF;
  else return st[i-1];
}
/*---------------------------------------------------------*/
double  readdouble(FILE *fpin)
{ 
  int j, *k; double f, *a; char *next;
  j=fscanf(fpin,"%lf",&f);
  if (j==1) return (f);
  else if (j==EOF) printout("error","readdouble: end of file, double not found\n");
  else  /* try to interpret as variable */
  { 
    next=readname(fpin);
    if (arraytype(next)=='d')
    { 
      if (arraysize(next)==1) {a=(double *)need(next); return (a[0]); }
    }
    else if (arraytype(next)=='i')
    { 
      if (arraysize(next)==1) {k=(int *)need(next); f=k[0]; return (f); }
    }
    printout("error","readdouble: input error:  next input item: %s\n",next);
    if (arraysize(next)>0)
      printout("error","readdouble expects a value or an array(type d or i) of size 1. size=%d, type=%c",
             arraysize(next),arraytype(next));
  }
  exitm4d(0);
  exit(1);  /* this statement not reached */
}
/*---------------------------------------------------------*/
char *readfilename(FILE *fpin)
{
  char *name, **namealias, *nameparse[3], intvalue[11];
  int i,jp[3]={0,0,0},n=0,*intvar;
  
  name=readnamep(fpin);
/* look for alias */
  if (arraytype(name)=='s')
  {
    namealias=(char **)need(name);
    printout("normal"," %s is alias for real file %s\n",name,namealias[0]);
    /* return (namealias[0]); */
    free(name);
    name=setname(namealias[0]);
  }
/* look to parse for integer variable */
  nameparse[0]=(char *)tmalloca(3*(strlen(name)+1),'c');
  nameparse[1]=nameparse[0]+strlen(name)+1;
  nameparse[2]=nameparse[1]+strlen(name)+1;
 
  for (i=0;i<strlen(name);i++)
  {
    if (name[i]=='#') n++;
    else 
    { 
      nameparse[n][jp[n]]=name[i];
      jp[n]++;
    }
    if (n>2) 
    {
      printout("warning readfilename","warning too many # in filename - no parsing\n");
      return (name);
    }
  }
  if (n==0) return (name); 
/* finish parsing */
  for (n=0;n<3;n++) nameparse[n][jp[n]]='\0';
  intvar=(int *)need(nameparse[1]);
  n=snprintf(intvalue,10,"%d",intvar[0]);
  intvalue[10]='\0';
  free(name);
  if (jp[2]==0) name=setname2(nameparse[0],intvalue);
  else name=setname3(nameparse[0],intvalue,nameparse[2]);
  return (name);
}
/*---------------------------------------------------------*/
float  readfloat(FILE *fpin)
{ 
  int j; float f; char *next;
  j=fscanf(fpin,"%f",&f);
  if (j==1) return (f);
  else 
  {
    printout("error","readfloat: float not found\n");
    if (j==EOF) printout("error readfloat"," end of file\n");
    else  
    { 
      next=readname(fpin);
      printout("error","  input error:  next input item: %s\n",next);
    } 
    exitm4d(0);
    exit(1); /* this statement not reached */
  }
}
/*---------------------------------------------------------*/
int  readint(FILE *fpin)
{
  int i,j, *k,n; char *next; double *a;
  j=fscanf(fpin,"%d",&i);
  if (j==1) return (i);
  else if (j==EOF) printout("error readint","readint: end of file, int not found\n");
  else  /* try to interpret as variable */
  { 
    next=readnamep(fpin);
    if (arraytype(next)=='i')
    { 
      if (arraysize(next)==1) {k=(int *)need(next); return (k[0]); }
    }
    if (arraytype(next)=='d')
    { 
      if (arraysize(next)==1) 
      {
        a=(double *)need(next); 
        n=a[0]+.01; if (a[0]<0) n=a[0]-.01;
        if (abs(n-a[0])<.01) return (n);
        printout("error reading","readint: error value of %s %lg not close enough to an integer\n",next,a[0]);
      }
    }
    printout("error readint","readint: input error:  next input item: %s\n",next);
    if (arraysize(next)>0)
      printout("error readinr","readint expects a value or an array(type d or i) of size 1. size=%d, type=%c",
             arraysize(next),arraytype(next));
  }
  exitm4d(0);
  exit(1); /* this statement not reached */
}
/*---------------------------------------------------------*/
char *readline(FILE *fpin)
{ 
  char c, namein[80], *name;
  int i, ii;
  i=0;
  while(1)
  { 
    c=getc(fpin);
    if (c==EOF && i==0) return (NULL);
    if (c != '\n' && c != '\r' && c != EOF && i<80) {  namein[i]=c; i++;}
    else 
    { 
      while (namein[i-1]==' ' && i>0) i--;  /* drop final blanks */
      name=(char *)tmalloca(i+1,'c');
      for (ii=0;ii<i;ii++) name[ii]=namein[ii];
      name[i]='\0';
      return (name);
    }
  }
}

/*---------------------------------------------------------*/   
char *readname(FILE *fpin)
{ 
  char c, cend, namein[202], *name;
  int i,ii;
  cend='a'; i=0;
  while(1)
  { 
    c=getc(fpin);
    if (c==EOF && i==0)  return (NULL);
    if (cend=='a') 
    { 
      if (c!=' ' && c!='\n' && c!='\t' && c!='\r' )
		{
        cend=' '; 
        if (c=='\'' || c=='\"') cend=c;
        if (cend==' ')
        {  
          namein[i]=c; 
          i++;
        }
		}
    }
    else if (c==cend || c=='\n' || c=='\r' || c=='\t' || c==EOF) 
    { 
      while (namein[i-1]==' ' && i>0) i--;  /* drop final blanks */ 
      name=(char *)tmalloca(i+1,'c');
      for (ii=0;ii<i;ii++) name[ii]=namein[ii];
      name[i]='\0';
      return (name);
    }
    else
    { 
      namein[i]=c; 
      i++;
      if (i==201)
      {
        namein[i]='\0';
        printout("error readname","error: name being read exceeds 200 characters\n %s\n",namein);
        exitm4d(0);
      }
    }
  }
}
/*---------------------------------------------------------*/   
char *readnamep(FILE *fpin)
{ 
  char c, cend, namein[202], *name;
  int i,ii;
  cend='a'; i=0;
  while(1)
  {
    c=getc(fpin);
    if (c==EOF && i==0)  return (NULL);
    if (cend=='a') 
    { 
      if (c!=' ' && c!='\n' && c!='\t' && c!='\r')
		{
        cend=' '; if (c=='\'' || c=='\"') cend=c;
        if (cend==' ')
        {  
          namein[i]=c; 
          i++;
        }
		}
    }
    else if (c==cend || c=='\n' || c=='\r' || c=='\t' || c==EOF) 
    { 
      while (namein[i-1]==' ' && i>0) i--;  /* drop final blanks */ 
      name=(char *)smalloc(sizeof(char)*(i+1));
      for (ii=0;ii<i;ii++) name[ii]=namein[ii];
      name[i]='\0';
      return (name);
    }
    else
    { 
      namein[i]=c; 
      i++;
      if (i==201)
      {
        namein[i]='\0';
        printout("error readnamep","error: name being read exceeds 200 characters\n %s\n",namein);
        exitm4d(0);
      }
    }
  }
}
/*-----------------------------------------------*/
char read1charname(FILE *fpin)
{ 
  char c, *name;
  name=readname(fpin);
  c=name[0];
  return c;
}

/*-----------------------------------------------*/
FILE *safefopen(const char *name, const char *how)
{
  FILE *fp;
  fp=fopen(name,how);
  if (fp==NULL) 
  {
    printout("error safeopen","error: cannot open file %s for %s\n",name,how);
    exitm4d(0);
  }
  return fp;
}

/* ---------------------------- */
char * setname(char*from)
{
  int i;
  char *name;
  i=strlen(from);
  name=(char *)smalloc(i+1);
  strcpy(name,from);
  return name;
}
/* ---------------------------- */
char * setname2(char*from,char*from2)
{
  int i;
  char *name;
  i=strlen(from)+strlen(from2);
  name=(char *)smalloc(i+1);
  strcpy(name,from);
  strcat(name,from2);
  return name;
}
/* ---------------------------- */
char * setname3(char*from,char*from2,char*from3)
{
  int i;
  char *name;
  i=strlen(from)+strlen(from2)+strlen(from3);
  name=(char *)smalloc(i+1);
  strcpy(name,from);
  strcat(name,from2);
  strcat(name,from3);
  return name;
}

