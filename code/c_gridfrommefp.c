/* contains c_gridfrommefp */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])
void readmefpgrid(char *filename, double **axyz, int **alt, double **aabc, int *idim);

void c_gridfrommefp(FILE *fpin, FILE *fprint)
{ 
  double *xyzmefp, *abcmefp; 
  int *ltmefp, idim[3];  /* mefp grid */
  char* filename;
  
  double *roundoff, tol=1.e-5;
  int *i4d; 
  double *xyz, *abcd;
  char *clt;
  
  int i,j,k,n,m,iall,imall,ii[3],nyreada[4]={1,1,1,1}; 
  double *aa,ff[3],flt,*abcm[3];
  char cclt;
  
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  filename=readfilename(fpin); /* read grid file name */
  readmefpgrid(filename,&xyzmefp,&ltmefp,&abcmefp,idim);  /* read mefp grid */
  abcm[0]=abcmefp; abcm[1]=abcm[0]+idim[0]; abcm[2]=abcm[1]+idim[1];
  
  i4d=(int*)createarray("idim4d",4,'i',0);
  Loop(k,0,4)   i4d[k]=readint(fpin); /* read grid dimensions */
  Loop(k,0,3)  /* if any of the first 3 dimensions are 0, copy from mefp-grid */
  {
    if (i4d[k]==0) { i4d[k]=idim[k], nyreada[k]=0; }
  }
  printout("normal","\ncalculation grid %d by %d by %d by %d\n",i4d[0],i4d[1],i4d[2],i4d[3]);
  
  abcd=(double*)createarray("abcd",i4d[0]+i4d[1]+i4d[2]+i4d[3],'d',0);
  aa=abcd;
  Loop(k,0,4)         /* read(or set) a,b,c,d for the grid */
  { 
    if (nyreada[k]==1)  Loop(i,0,i4d[k]) aa[i]=readdouble(fpin);  
    else  Loop(i,0,i4d[k]) aa[i]=abcm[k][i];
    printout("normal","dim %d parm=  %lg",i4d[k],aa[0]);
    Loop(i,1,i4d[k])
    { 
      if (i%10==0) printout("normal","\n          ");
      printout("normal"," %lg",aa[i]);
      if (aa[i]<aa[i-1]) printout("warning c_gridfrommefp","  warning, point out of order\n");
    }
    printout("normal","\n");
    aa+=i4d[k];
  }
  
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  imall=idim[0]*idim[1]*idim[2];
  xyz=(double*)createarray("xyz",iall*3,'d',0);
  clt=(char*)createarray("clt",iall,'c',0);
  
  for (i=0;i<i4d[0];i++) for (j=0;j<i4d[1];j++) for (k=0;k<i4d[2];k++)
  {  
    m=i+i4d[0]*(j+i4d[1]*(k+i4d[2]));
    ii[0]=findex(abcd[i],abcmefp,idim[0],ff);
    ii[1]=findex(*(abcd+i4d[0]+j),abcmefp+idim[0],idim[1],ff+1);
    ii[2]=findex(*(abcd+i4d[0]+i4d[1]+k),abcmefp+idim[0]+idim[1],idim[2],ff+2);
    flt=Value3d(ltmefp,idim,ii,ff);
    cclt='w';   	   if (flt<3-tol) cclt='f';   if (flt>3+tol) cclt='s';
    for (n=0;n<i4d[3];n++) 
    {  
      m=i+i4d[0]*(j+i4d[1]*(k+i4d[2]*n));
      xyz[m]=Value3d(xyzmefp,idim,ii,ff);
      xyz[m+iall]=Value3d((xyzmefp+imall),idim,ii,ff);
      xyz[m+2*iall]=Value3d((xyzmefp+2*imall),idim,ii,ff);
      /*  xyz[m+3*iall]=abcd[i4d[0]+i4d[1]+i4d[2]+n];  for old xyzt array */
      clt[m]=cclt;
    }
  }
}
