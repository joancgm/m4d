/* contains c_ppts2eqppts */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_rmsminmax(FILE *fpin, FILE *fprint)
{ 
  int *i4d;  /* need */
  double *v; int *ITER; double *TIME; /* use if available */
  char *name; /* input parm */
  int i4dm[4],i4dp[4],i4ddum[4]={1,1,1,1},iallg,iallm,iallp;
  int *ideq,size,i,ii[4],imin,imax;
  double vmin,vmax,vrms,*rmsmm;
  int n=100,isave=10,k=0;
  
  i4d=(int *)need("idim4d");
  ITER=(int *)find("ITER");
  TIME=(double *)find("TIME");
  rmsmm=(double *)createarray("rmsmm",2*isave,'d',0);
  Loop(i,0,2*isave) rmsmm[i]=0;
  
  Loop(i,0,3) {i4dp[i]=i4d[i]+1; i4dm[i]=i4d[i]-1;} i4dp[3]=i4d[3]; i4dm[3]=i4d[3];
  iallg=Prod4(i4d); iallp=Prod4(i4dp); iallm=Prod4(i4dm);
  
  while (n)
  { 
    name=readname(fpin);
    if (name[0]=='\0')  break;
    if (strcmp(name,"c:")==0) { fseek(fpin,(long)(-3),1); break; }
    n--;	
    size=arraysize(name);
    if (size==0) {printout("normal"," array %s not found\n",name); continue;}
    if (size%iallg==0) ideq=i4d;
    else if (size%iallp==0) ideq=i4dp;
    else if (size%iallm==0) ideq=i4dm;
    else {i4ddum[0]=size; ideq=i4ddum; }
    
    v=(double *)need(name);
    vmin=v[0]; vmax=v[0]; imin=0; imax=0;
    vrms=0;
    Loop(i,0,size)
    { 
      vrms+=v[i]*v[i];
      if (v[i]<vmin) {vmin=v[i]; imin=i; }
      if (v[i]>vmax) {vmax=v[i]; imax=i;}
    }
    vrms=sqrt(vrms/size);
    printout("normal","RMSMM %s rms %lg,",name,vrms);
    iexpand(imin,ideq,ii);
    printout("normal"," min %lg at %d %d %d %d, ",vmin,ii[0],ii[1],ii[2],ii[3]);
    iexpand(imax,ideq,ii);
    printout("normal"," max %lg at %d %d %d %d, ",vmax,ii[0],ii[1],ii[2],ii[3]);		
    if (ITER!=0) printout("normal"," ITER %d ",*ITER);
    if (TIME!=0) printout("normal"," TIME %lg",*TIME);
    printout("normal","\n");
    /* save first isave to rmsmm */
    if (k<isave) { rmsmm[2*k]=vrms; rmsmm[2*k+1]=max(vmax,-vmin); k++;	}
  }
}


