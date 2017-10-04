/* contains c_omrhscoakley */
#include "global.h"
/* rhs of omturb equation, Coakley model */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_omrhscoakley(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*noindfwpts,*wherep,*wherefw; /* needed arrays */
  double *q, *om, *rho, *vlam, *pkdk, *walldist;
  char *cltp;
  int *itrange; /* use if available */
  double *rhsom, *caom;  /* set rhsom */
  
  int jcoef=0; double fac=1; char oldnew='n'; /* fixed parm */
  
  int i,k,its,ite,ip[4],i4dp[4],i4dm[4],ii[4],mc[8],ia[3],mmid,iallm,mwp;
  double cave[8],vol[8],rhom,qm,omm,p,d,ddo,ce1,ce2,fv,yqdnu;
  int two[3]={2,2,2};
  
  /* jcoef=readint(fpin); omit as parameter, uncomment this to reinstate */
  printout("normal","coefficient set %d used for convection, viscous, and/or time\n",jcoef);
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  noindfwpts=(int *)need("noindfwpts");
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  wherefw=(int *)need("wherefw");
  q=(double *)need("qturb"); 
  om=(double *)need("omturb"); 
  rho=(double *)need("rho");
  walldist=(double *)need("walldist");
  vlam=(double *)need("vlam");
  pkdk=(double *)need("pkdk");
  
  if (jcoef>=0)    /* convection etc */
    coefrhs(fprint,"omturb","rhsomturb",jcoef,oldnew,fac);
  
  rhsom=(double *)need("rhsomturb");
  caom=(double *)createarray("caomturb",noindfwpts[0],'d',0);
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2]) caom[i]=0;
  
  geom8init();
  geomcinit();
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* each cont. c.v. */
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
    mwp=wherep[In4(ip,i4dp)];
    if (mwp<0) continue;  /* not valid c.v. leave as zero */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
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
    
    yqdnu=walldist[mwp]*qm/(vlam[mmid]/rhom);
    fv=.02*yqdnu;
    if (fv>1.e-5) fv=1-exp(-fv);
    ce1=1.555;
    if (fv>0) ce1=1.5+.055/fv;
    ce2=1.833;
    p=(ce1-1)*pkdk[mmid]*rhom;
    ddo=(ce2-1)*rhom;
    Loop(k,0,8)
    {
      d=ddo*max(.5*om[mc[k]],min(2.*om[mc[k]],omm));
      if (p>0) rhsom[wherefw[mc[k]]]+= vol[k]*(p*omm-d*om[mc[k]]);
      else rhsom[wherefw[mc[k]]]+= vol[k]*(p*om[mc[k]]-d*om[mc[k]]);
      caom[wherefw[mc[k]]]+=vol[k]*(abs(p)+2*abs(d));
    }
  } /* each cont c.v. */
}	 
