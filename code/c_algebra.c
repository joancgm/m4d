/* contains c_algebra */
#include "global.h"

void c_algebra(FILE *fpin, FILE *fprint)

/*  read namea nameb namec named e f g
 a=b(e*c+f)**g+d
 */
{ 
  char *name[4];
  int i,imin,imax, n;
  double e,f,g,vv, *vout,*vp[4];;
  Array *v[4];
  
  imin=0; imax=0;
  for (n=0;n<4;n++) 
  { 
    name[n]=readname(fpin); 
    v[n]=(Array *)findarray(name[n]);
    /* check for array, exit on not found unless empty name */
    if ( n>0 && v[n]==0 && name[n][0]!='\0')
    { 
      printout("error c_algebra"," array %s not found, exit\n",name[n]); 
      exitm4d(0); 
    }
    if (v[n]!=0)
    { 
      if (imin==0) {imin=v[n]->size; imax=imin;}
      else {imin=min(imin,v[n]->size); imax=max(imax,v[n]->size);}
	   if (v[n]->type != 'd') 
      {
        printout("error c_algebra"," array %s not of type 'd'\n",name[n]); 
        exitm4d(0); 
      }
	   vp[n]=(double *)v[n]->pointer;
    }
  }
  e=readdouble(fpin);
  f=readdouble(fpin);
  g=readdouble(fpin);
  printout("normal","%s = %s *(%lg * %s + %lg)^%lg + %s\n",
          name[0],name[1],e,name[2],f,g,name[3]);
  if (imax>imin) printout("warning c_algebra","warning arrays have different sizes %d %d\n",imin,imax);
  if (v[0]==0) 
  { 
    if (imin>0)  vout=(double *)createarray(name[0],imin,'d',0);
    else 
    { 
      printout("warning c_algebra"," no arrays found, algebra ignored\n");
      return;
    }
  }
  else vout=(double *)v[0]->pointer;
  
  for (i=0;i<imin;i++)
  { 
    vv=1.;
    if (v[2]!=0) vv *= vp[2][i];
    vv = e * vv + f;
    if (g !=1. && vv>0) vv=pow(vv,g);
    if (v[1]!=0) vv *= vp[1][i];
    if (v[3]!=0) vv += vp[3][i];
    vout[i] = vv;
  }
}


