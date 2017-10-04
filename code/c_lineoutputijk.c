/* contains c_lineoutputijk */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/*  lineoutput of dimensioned variables can do diagonal as well as flat grid line output */
void c_lineoutputijk(FILE *fpin, FILE *fprint)
{ 
  char *filename, *name[20], nameif[10];   /* input */
  int nf, idf[10], isf[10], ief[10]; 
  int *i4d;    /* use if available */
  
  int n, i,k, ifacf[10], sign[10], ivalues, numvar=0, ii[10], ipt;
  double *v[20];
  FILE *fline;
  
  i4d=(int*)find("idim4d");
  filename=readfilename(fpin);                  /* output file */
  printout("normal","lineoutputijk to file %s:",filename);
  fline=safefopen(filename,"w");  
  
  Loop(i,0,20)           /* variables */
  {
    name[i]=readname(fpin);
    if (name[i][0]=='\0') break; 
    printout("normal"," %s",name[i]);
    v[i]=(double*)find(name[i]);
    if(v[i]==0) {i--; printout("normal"," not found,"); }
    numvar=i+1;
  }
  printout("normal","\n");
  
  nf=readint(fpin);    /* dimensions, same for all variables! */
  printout("normal"," dims %d\n",nf);
  if (nf>10) {printout("warning c_lineoutoutijk"," routine only dimensioned for 10 dimensions\n"); return; }
  for (n=0;n<nf;n++)
  { 
    nameif[n]=read1charname(fpin);
    idf[n]=readint(fpin);
    if (i4d>0)
    {
      if (nameif[n]=='i' ) idf[n] += i4d[0];
      else if (nameif[n]=='j') idf[n] += i4d[1];
      else if (nameif[n]=='k') idf[n] += i4d[2];
      else if (nameif[n]=='t') idf[n] += i4d[3];
    }
    printout("normal"," %c %d",nameif[n],idf[n]);
    if (n==0) ifacf[n]=1; else ifacf[n]=ifacf[n-1]*idf[n-1];
  }
  /* output file header */
  fprintf(fline,"m4d lineoutputijk\n%d\n ",nf+numvar);
  Loop(n,0,nf) fprintf(fline," %c",nameif[n]); 	fprintf(fline,"\n");
  Loop(i,0,numvar) fprintf(fline," %s",name[i]); 	fprintf(fline,"\n");
  
  while(1) /* each range to be output set first value to 0 to end */
  {
    printout("normal"," range");
    ivalues=0;
    for (n=0;n<nf;n++)
    {
      isf[n]=readint(fpin);
      isf[n+1]=isf[n];
      if (isf[n]<=0) break;
      ief[n]=readint(fpin);
      isf[n]=min(idf[n],isf[n]);
      ief[n]=min(idf[n],max(1,ief[n]));	
      printout("normal"," %c %d %d",nameif[n],isf[n],ief[n]);
      sign[n]=1; if (ief[n]<isf[n]) sign[n]=-1;
      ivalues=max(sign[n]*(ief[n]-isf[n])+1,ivalues);
    }
    if (isf[n]<=0) break;
    printout("normal",", %d values\n",ivalues);
    
    Loop(n,0,nf) ii[n]=isf[n]-1;
    Loop(k,0,ivalues)
    {
      if (k==0) Loop(n,0,nf) ii[n]=isf[n]-1;
      else Loop(n,0,nf) { ii[n]+=sign[n]; if (ii[n]==ief[n]) ii[n]=isf[n]-1; }
      ipt=ii[0]; Loop(n,1,nf) ipt+=ii[n]*ifacf[n];
      Loop(n,0,nf) fprintf(fline," %d",ii[n]);
      Loop(i,0,numvar)  fprintf(fline," %lg",v[i][ipt]);
      fprintf(fline,"\n");
    }
  }
  fclose(fline);
}
