/* contains c_lineoutput */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* read filename for line output in lineplot data format
 names of variables ending with "" (max 20) on the points variables or line variables
 clinedir: one of  a b c d  i j k n 
 starting location (4 values) a,b,c,d or i,j,k,n
 end location 1 value
 repeat from read clinedir until clinedir not an allowable value (suggest e for end)
 */
void c_lineoutput(FILE *fpin, FILE *fprint)
{ 
  char *filename, *name[20], *clinedir;   /* input */
  double ast[4],abcend;
  int ist[4],ijkend;
  
  int *i4d;    /* need */
  double *a[4],*v[20];
  
  FILE *fline;
  int i,j,numvar=0,nya=0,idir=0,di,ii[4],ipt,iall;
  double f[1];
  char ijkn[5]="ijkn";
  
  i4d=(int*)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double*)need("abcd");
  Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
  
  filename=readfilename(fpin);                  /* input */
  printout("normal","line output to file %s: for ",filename);
  fline=safefopen(filename,"w");           /* write file */
  Loop(i,0,20)
  { 
    name[i]=readname(fpin);
    if (name[i][0]=='\0') break; 
    printout("normal"," %s",name[i]);
    v[i]=(double*)find(name[i]);
    if(v[i]==0) {i--; printout("warning c_lineoutput"," not found,"); }
    else if (arraysize(name[i])%iall != 0) {i--; printout("warning c_lineoutput"," incompatible size,"); } 
    numvar=i+1;
  }
  printout("normal","\n");
  
  fprintf(fline,"m4d lineoutput\n %d\n abcd idir",numvar+2);
  Loop(i,0,numvar) fprintf(fline," %s",name[i]); 	fprintf(fline,"\n");
  
  while(1)
  { 
    clinedir=readname(fpin);
    if (strcmp(clinedir,"c:")==0){fseek(fpin,(long)(-3),1); break; }
    if (clinedir[0]=='a') {nya=1; idir=0;}
    else if (clinedir[0]=='b') {nya=1; idir=1;}
    else if (clinedir[0]=='c') {nya=1; idir=2;}
    else if (clinedir[0]=='d') {nya=1; idir=3;}
    else if (clinedir[0]=='i') {nya=0; idir=0;}
    else if (clinedir[0]=='j') {nya=0; idir=1;}
    else if (clinedir[0]=='k') {nya=0; idir=2;}
    else if (clinedir[0]=='n') {nya=0; idir=3;}
    else break;
    if (nya==1)
    { 
      Loop(i,0,4) ast[i]=readdouble(fpin);
      abcend=readdouble(fpin);
      printout("normal","line from abcd=%lg,%lg,%lg,%lg to %c=%lg",
              ast[0],ast[1],ast[2],ast[3],clinedir[0],abcend);
      Loop(i,0,4) ist[i]=findex(ast[i],a[i],i4d[i],f)+f[0]+.5;
      ijkend=findex(abcend,a[idir],i4d[idir],f)+f[0]+.5;
      printout("normal"," (line from ijkn=%d %d %d %d to %c=%d)\n",
              ist[0],ist[1],ist[2],ist[3],ijkn[idir],ijkend);
    }
    else
    { 
      Loop(i,0,4) ist[i]=readint(fpin);
      ijkend=readint(fpin);
      if (ist[0]>0 || ist[1]>0 || ist[2]>0 || ist[3]>0) /* limit if really 4-d arrray */
      { 
        Loop(i,0,4) ist[i]=min(ist[i],i4d[i]);
        ijkend=min(ijkend,i4d[idir]);
      }
      printout("normal","line from ijkn=%d %d %d %d to %c=%d\n",
              ist[0],ist[1],ist[2],ist[3],ijkn[idir],ijkend);
    }
    di=1; if (ijkend<ist[idir]) di=-1;
    Loop(i,0,4) ii[i]=ist[i];
    Loop(i,0,abs(ijkend-ist[idir])+1)
    { 
      ipt=In4(ii,i4d);
      if (nya==1) fprintf(fline,"%lg %d",a[idir][ii[idir]],idir);
      else fprintf(fline,"%d %d",ii[idir],idir);
      Loop(j,0,numvar) fprintf(fline," %lg",(double)v[j][ipt]);
      fprintf(fline,"\n");
      ii[idir]+=di;
    }
  }
  fclose(fline);
}
