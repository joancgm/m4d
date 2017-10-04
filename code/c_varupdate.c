/* contains c_varupdate */
#include "global.h"
#include "fixbij.h"

/* update a variable with the option of limiters and method 
 type= l log on neg changes, 
 tyoe= i inverse on negative changes
 type= b bij update linearly then apply realizability
 
 */
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_varupdate(FILE *fpin, FILE *fprint)
{ 
  int *i4d;  /* need */
  double *v, *dv;
  int *itrange; /* use if available*/
  char *name,*dname,*type; /* input parameters */
  double *updatew; /* create */
  
  int i,j,iall,iall3,ist,iend,irep,iwarn=0,iwarnz=0,iwarnh=0; 
  double fac,facmin=1;
  double bij[6],dbij[6],bijn[6],ch,vnew,dbija;
  double gtol=.001;
  int bug=10;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d); 
  iall3=i4d[0]*i4d[1]*i4d[2];
  itrange=(int *)find("itrange");
  if (itrange>0) { ist=itrange[0]*iall3; iend=(itrange[1]+1)*iall3;}
  else {ist=0; iend=iall; }
  
  name=readname(fpin);
  dname=readname(fpin);
  type=readname(fpin);
  v=(double *)need(name);
  dv=(double *)need(dname);
  printout("normal","update %s with %s method %s\n",
          name,dname,type);
  irep=arraysize(name)/iall;
  updatew=(double *)createarray("updatew",3,'d',0);
  
  if (type[0]=='l' || type[0]=='i') 
    Loop(i,ist,iend)
  {
    if (dv[i]>=0) v[i]+=dv[i];
    else if (v[i]>0) 
    {
      if (type[0]=='l') vnew=v[i]*exp(dv[i]/v[i]);
      else vnew=v[i]*1./(1.-dv[i]/v[i]);
      ch=vnew-v[i];
      v[i]=vnew;
      fac=ch/dv[i];
      facmin=min(fac,facmin);
      if (fac<.5) iwarnh++;
      if (fac<.001) iwarnz++;
    }
    else {v[i]=0; iwarnz++;}
  }
  
  else if (type[0]=='b' && irep==6)  /* v=bij */
    Loop(i,ist,iend)
  {
    Loop(j,0,6) 
    {
      bij[j]=v[i+j*iall]; 
      dbij[j]=dv[i+j*iall]; 
      bijn[j]=bij[j]+dbij[j];
    }
    fac=fixbij(bijn,gtol);
    if (fac<1) iwarn++;
    Loop(j,0,6) v[i+j*iall]=bijn[j];
    if (fac<.9 && bug>0)
    { 
      printout("normal","%d  bij %lg %lg %lg %lg %lg %lg\n",
              i,bij[0],bij[1],bij[2],bij[3],bij[4],bij[5]);
      printout("normal","%d dbij %lg %lg %lg %lg %lg %lg\n",
              i,dbij[0],dbij[1],dbij[2],dbij[3],dbij[4],dbij[5]);
      printout("normal","%d bijn %lg %lg %lg %lg %lg %lg\n",
              i,bijn[0],bijn[1],bijn[2],bijn[3],bijn[4],bijn[5]);
      bug--;
    } 
    
    dbija=0; ch=0;
    Loop(j,0,6)
    {
      dbija+=abs(dbij[j]);
      ch+=abs(bij[j]-bijn[j]);
    }
    fac=1;
    if (dbija>0)
    {
      fac=ch/dbija;
      if (fac<.5) iwarnh++;
      if (fac<.001) iwarnz++;
    }
    facmin=min(fac,facmin);
  }
  else 
  {
    printout("error c_varupdate","Error could not update method %c, irep %d\n",type[0],irep);
    exitm4d(0);
  }
  updatew[0]=facmin;
  updatew[1]=iwarnh;
  updatew[2]=iwarnz;
  printout("normal","update %s: facmin %lg, fac<.5 for %d pts, <.001 for %d pts\n",name,facmin,iwarnh,iwarnz);
  
}