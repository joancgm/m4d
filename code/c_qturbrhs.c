/* contains c_qturbrhs */
#include "global.h"
/* rhs of qturb equation, generic using P/k  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_qturbrhs(FILE *fpin, FILE *fprint)
{
  int *i4d,*noindfwpts,*wherep,*wherefw; /* needed arrays */
  double *q, *om, *rho, *pkdk;
  char *cltp;
  int *itrange; /* use if available */
  double *rhsq, *caq;  /* set rhsq */
  
  int jcoef=0; double fac=1; char oldnew='n'; /*fixed parm */
  
  int i,k,its,ite,ip[4],i4dp[4],i4dm[4],ii[4],mc[8],ia[3],mmid;
  double cave[8],vol[8],rhom,qm,omm,p,d;
  int two[3]={2,2,2};
  
  /* jcoef=readint(fpin); omit as parameter, uncomment this to reinstate */
  printout("normal","coefficient set %d used for convection, viscous, and/or time\n",jcoef);
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  noindfwpts=(int *)need("noindfwpts");
  wherep=(int *)need("wherep");
  wherefw=(int *)need("wherefw");
  cltp=(char *)need("cltp");
  q=(double *)need("qturb"); 
  om=(double *)need("omturb"); 
  rho=(double *)need("rho");
  pkdk=(double *)need("pkdk");
  if (jcoef>=0)    /* convection etc */
    coefrhs(fprint,"qturb","rhsqturb",jcoef,oldnew,fac);
  
  rhsq=(double *)need("rhsqturb");
  caq=(double *)find("caqturb");
  if (caq==0) 
  { 
    caq=(double *)createarray("caqturb",noindfwpts[0],'d',0);
    Loop(i,0,noindfwpts[0]) caq[i]=0;
  }
  else Loop(i,noindfwpts[its+1],noindfwpts[ite+2]) caq[i]=0;
  
  geom8init();
  geomcinit();
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* each cont. c.v. */
  {
    if (wherep[In4(ip,i4dp)]<0) continue; /* do only for valid volumes */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
    geomcvave(cave,ip);
    geom8vol(vol,ip);
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    {
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }
    omm=0; qm=0; rhom=0; /* average values for rho q and om */
    Loop(k,0,8)  
    { 
      rhom+=cave[k]*rho[mc[k]];
      /* omm+=cave[k]/om[mc[k]]; */ /* assume 1/om varies linearly */
      omm+=cave[k]*log(om[mc[k]]);
      qm+=cave[k]*q[mc[k]];
    }
    /*omm=1./omm; */
    omm=exp(omm);     
    
    p=0.5*pkdk[mmid]*rhom;
    d=.5*omm*rhom;
    Loop(k,0,8)
    { 
      if (p>0) rhsq[wherefw[mc[k]]]+= vol[k]*(p*qm-d*q[mc[k]]);
		else rhsq[wherefw[mc[k]]]+= vol[k]*(p*q[mc[k]]-d*q[mc[k]]);
      caq[wherefw[mc[k]]]+=vol[k]*(abs(p)+abs(d));
    }
  } /* each cont c.v. */
}	 
