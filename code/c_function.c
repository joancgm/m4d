/* contains c_function */
#include "global.h"

void c_function(FILE *fpin, FILE *fprint)

/*  functions: v0=arctan2(v1,v2)  artan2deg,  v0=tan(v1)  tandeg,
 v0=log(v1) log10  v0=exp(v1)
 v0=max(v1,v2) v0=min(v1,v2)  v0=abs(v1) sin(v1);  cos(v1);
 tanh(v1)
 */
{ 
  char *name[4], *namefunction;
  
  int  i,j,n,size[4],sizemax,nmax,nneg,nzero,nynew=0;
  double *v[4]; 
  double pi=3.141592654;
  
  name[0]=readname(fpin); 
  namefunction=readname(fpin);
  name[1]=readname(fpin);
  nmax=2;
  /*       extra name for arctan2, max, min */
  if (strcmp(namefunction,"arctan2")==0 
      || strcmp(namefunction,"arctan2deg")==0
      || strcmp(namefunction,"max")==0
      || strcmp(namefunction,"min")==0)
  { 
    nmax=3;
    name[2]=readname(fpin); 
  } 
  printout("normal","%s=%s(%s",name[0],namefunction,name[1]);
  if (nmax>2) printout("normal",",%s",name[2]);
  printout("normal",")\n");
  
  /* find needed arrays */
  sizemax=0;
  for (n=1;n<nmax;n++) 
  { 
    v[n]=(double *)need(name[n]);
    size[n]=arraysize(name[n]);
    if (size[n]>sizemax) sizemax=size[n];
  }
  /* find or create target array */
  v[0]=(double *)find(name[0]);
  if (v[0]==0) 
  { 
    v[0]=(double *)createarray(name[0],sizemax,'d',0);
    printout("normal","array %s created with size %d\n",name[0],sizemax);
    nynew=1;
  }
  size[0]=arraysize(name[0]);
  /* check sizes */
  for (n=1;n<nmax;n++) 
    if (size[n]!=size[0] && size[n]!=1)
      printout("normal","size warning %s size %d,  %s size %d\n",
              name[0],size[0],name[n],size[n]);
  
  if (strcmp(namefunction,"arctan2")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=atan2(v[1][min(i,size[1]-1)],v[2][min(i,size[2]-1)]);
  
  else if (strcmp(namefunction,"arctan2deg")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=atan2(v[1][min(i,size[1]-1)],v[2][min(i,size[2]-1)])*180./pi;
  
  else if (strcmp(namefunction,"tan")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=tan(v[1][min(i,size[1]-1)]);
  
  else if (strcmp(namefunction,"sin")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=sin(v[1][min(i,size[1]-1)]);
  
  else if (strcmp(namefunction,"cos")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=cos(v[1][min(i,size[1]-1)]);
  
  else if (strcmp(namefunction,"tanh")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=tanh(v[1][min(i,size[1]-1)]);
  
  else if (strcmp(namefunction,"tandeg")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=tan(v[1][min(i,size[1]-1)])*pi/180.;
  
  else if (strcmp(namefunction,"log")==0)
  {  
    nneg=0; nzero=0;
    for (i=0;i<size[0];i++) 
    { 
      j=min(i,size[1]-1);
      if (v[1][j]>0) v[0][i]=log(v[1][j]);
      else if (v[1][j]<0) { nneg++; v[0][i]=log(-v[1][j]);}
      else if (v[1][j]==0) { nzero++; v[0][i]=-50.;}
    }
    if (nneg>0) printout("warning c_function","log(%s<0) set to log(-%s) at %d points\n",
                        name[1],name[1],nneg);
    if (nzero>0) printout("warning c_function","log(%s=0) set to -50 at %d points\n",
                         name[1],nzero);
  }
  else if (strcmp(namefunction,"log10")==0)
  { 
    nneg=0; nzero=0;
    for (i=0;i<size[0];i++) 
    { 
      j=min(i,size[1]-1);
      if (v[1][j]>0) v[0][i]=log10(v[1][j]);
      else if (v[1][j]<0) { nneg++; v[0][i]=log10(-v[1][j]);}
      else if (v[1][j]==0) { nzero++; v[0][i]=-50.;}
    }
    if (nneg>0) printout("warning c_function","log10(%s<0) set to log10(-%s) at %d points\n",
                        name[1],name[1],nneg);
    if (nzero>0) printout("warning c_function","log10(%s=0) set to -50 at %d points\n",
                         name[1],nzero);
  }
  
  else if (strcmp(namefunction,"exp")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=exp(v[1][min(i,size[1]-1)]);
  
  else if (strcmp(namefunction,"max")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=max(v[1][min(i,size[1]-1)],v[2][min(i,size[2]-1)]);
  
  else if (strcmp(namefunction,"min")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=min(v[1][min(i,size[1]-1)],v[2][min(i,size[2]-1)]);
  
  else if (strcmp(namefunction,"abs")==0)
    for (i=0;i<size[0];i++) 
      v[0][i]=abs(v[1][min(i,size[1]-1)]);
  
  else
  {
    printout("warning c_function","function %s not implemented, \n  function=arctan2 arctan2deg tan tandeg log log10 max min abs\n",namefunction);
    if (nynew==1) 
    {
      for (i=0;i<size[0];i++) v[0][i]=0.;
      printout("warning"," %s set to zero\n",name[0]);
    }
  }
}

