/* contains m4dtitle main printout exitm4d */
#include "global.h"
#include <sys/times.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "commandlist.h"
#ifndef RUSAGE_SELF
#define RUSAGE_SELF 0
#endif

typedef struct { char *name; FILE *fpin; int loop; } Infile;

void m4dtitle(FILE *fprint)
{ 
  printout("normal","\nM4D - a research CFD code by Joan G. Moore\n");
  printout("normal","  for steady or time accurate calculations\n");
  printout("normal","  for inviscid, laminar or turbulent flow\n");
  printout("normal","              Features: \n");
  printout("normal","    Convection adapted control volumes\n");
  printout("normal","The transitional MARVS Reynolds stress model\n");
  printout("normal","            Version 2015.2c\n\n");
}
/* -----------------------------  main -------------------------*/
int  main(void)
{
  FILE *fpin,*fprint;
  char *command,c;
  char namestdin[]="stdin";
  char namestdout[]="stdout";
  char nameon[]="on";
  char nameoff[]="off";
  Infile infiles[20]; int innow; /* infile information */
  struct tms time; double tstart,tnow,tlast; /* cpu time tracking */
  struct rusage rus; long rumaxrss;  /* max space tracking */
  double *timespace;   /* put into dumpable array */
  char **info4print;   /* save stuff for print management */
  
  innow=0;
  infiles[0].fpin=stdin; infiles[0].loop=1; infiles[0].name=namestdin;
  
  fpin=infiles[innow].fpin;
  fprint=stdout;
  times(&time);
  tstart=.01*time.tms_utime;
  tnow=0;
  tlast=0;
  info4print=(char **)createarray("info4print",20,'p',0);
  info4print[1]=namestdin;
  info4print[2]=(char *)fprint;
  info4print[3]=namestdout;
  info4print[4]=nameon;
  timespace=(double *)createarray("timespace",2,'d',0);
  
  m4dtitle(fprint);    /* print title */
  
  while (1)
  {
    c=findstring(fpin,"c:"); 
    if (c==EOF)    /* stacked input files */
    { 
      if (infiles[innow].loop>1)
		{ rewind(fpin); infiles[innow].loop--; }
		else 
		{ 
        if (innow==0){printout("normal","end of file\n"); exitm4d(0);}
        fclose(fpin);
        free(infiles[innow].name);
        innow--;
        fpin=infiles[innow].fpin;
		}
      printout("normal","\n--------------input from file %s, loops remaining %d\n\n",
               infiles[innow].name,infiles[innow].loop);
      continue;
    }
    command=readname(fpin);
    info4print[0]=command;
    times(&time); tnow=.01*time.tms_utime-tstart;
    getrusage(RUSAGE_SELF,&rus);
    rumaxrss=rus.ru_maxrss;
    timespace[0]=tnow;
    timespace[1]=rumaxrss;
    
    printout("normal","\n----- c: %s -----     time %g   change %g  maxsize %ld\n",
             command,tnow,tnow-tlast,rumaxrss);
    tlast=tnow;
    
    /* -------------------command end -----------------------*/
    if (strcmp(command,"end")==0) exitm4d(0);
    
    /* --------------- command infile ---------------------*/
    else if (strcmp(command,"infile")==0)  /* add infile to stack and read it */
    {
      innow++; 
      if (innow>19) {printout("error infile (main)","error, 20 stacked infiles max\n"); exitm4d(0); }
      infiles[innow].name=readfilename(fpin);  /* read file name */
      printout("normal","\n--------------input from file %s", infiles[innow].name);
      infiles[innow].loop=readint(fpin);   /* read loops */
      if (infiles[innow].loop<1) infiles[innow].loop=1;
      printout("normal",", loops remaining %d\n\n",infiles[innow].loop);
      infiles[innow].fpin=safefopen(infiles[innow].name,"r");
      fpin=infiles[innow].fpin;
    }
    /* ----------------------command printcontrol -------------------------*/
    /*                   works with subroutine printout below          */
    else if (strcmp(command,"printcontrol")==0)  /* print instructions or files */
    {
      char *action,*newfile; FILE *newfprint;
      action=readname(fpin);
      if (strcmp(action,"on")==0) 
      {
        if (info4print[4]==nameoff) 
          fprintf(fprint,"\n printcontrol: resuming print on %s\n",info4print[3]);
        info4print[4]=nameon;
      }
      else if (strcmp(action,"off")==0)
      {
        if (info4print[4]==nameon)
          fprintf(fprint,"\n printcontrol: stopping normal print on %s\n",info4print[3]);  
        info4print[4]=nameoff;
      }
      else if ((strcmp(action,"filenew")==0) || (strcmp(action,"fileold") ==0))
      {
        newfile=readfilename(fpin);
        fprintf(fprint,"\n printcontrol: changing print from %s to %s\n",info4print[3],newfile);
        if (strcmp(newfile,namestdout)==0)
        {
          if (info4print[3]!=namestdout) { fclose((FILE *)info4print[2]); free(info4print[3]); }
          info4print[2]=(char *)stdout;
          info4print[3]=namestdout;
        }
        else /* real file */
        {
          if (strcmp(action,"filenew")==0) newfprint=safefopen(newfile,"w");
          else newfprint=safefopen(newfile,"a");
          if (info4print[3]!=namestdout) { fclose((FILE *)info4print[2]); free(info4print[3]); }
          info4print[2]=(char *)newfprint;
          info4print[3]=newfile;
          fprintf(newfprint,"\n printcontrol: m4d print file ...  \n");
        }
        fprint=(FILE *)info4print[2];
      }
      else printout("warning printcontrol (main)"," printcontrol action %s not valid, ignored\n",action);
    }
    
    /* ------ the remaining subroutine commands are listed in commandlist.c ------ */
    else if (commandlist(command,fpin,fprint)==0) 
    { printout("error main"," command %s not found, exit\n",command); exitm4d(0); }
    tmalloca(-1,'c');
  }
}

/*     -------------------------- subroutine printout ------------------*/
/* info4print: 0-command, 1-infilename, 2-print file pointer, 3-printfilename
 4-status (on or off)
 */
int printout(const char *type, const char *format, ...)
{
  va_list argp;
  char **info4print;
  int i=0;
  
  info4print=(char **)need("info4print");
  if (strcmp(info4print[4],"on") ==0)
  {
    va_start(argp, format);
    i=vfprintf((FILE *)info4print[2],format,argp);
    va_end(argp);
    return (i);
  }
  /* normal print is off */
  if (type[0]=='n') return 0;
  if ((type[0]=='w' && strlen(type)>7) || 
      (type[0]=='e' && strlen(type)>5) )
  {
    fprintf((FILE *)info4print[2],"%s, command=%s, infile=%s\n",
            type,info4print[0],info4print[1]);
  }
  va_start(argp, format);
  i=vfprintf((FILE *)info4print[2],format,argp);
  va_end(argp);
  return (i);
}
/*   ------------------------- exitm4d ----------------------------*/
void exitm4d(int i)
{
  /* putting exit through here in case want to do something */
  exit(0);
}
