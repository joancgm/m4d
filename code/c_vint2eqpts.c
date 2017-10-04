/* contains c_vint2eqpts */
#include "global.h"
/* between pts name, volume integrated to eq pts namevint  (i.e. contribution to r.h.s. of eqn) */
/* read name namevint oldnew(char) */
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_vint2eqpts(FILE *fpin, FILE *fprint)
{
  int *i4d,*noindfwpts,*wherep,*wherefw; /* needed arrays */
  double *v;
  char *cltp;
  int *itrange; /* use if available */
  double *vint;  /* set or update */
  
  char *name, *namevint, oldnew;
  
  int i,k,its,ite,ip[4],i4dp[4],i4dm[4],ii[4],mc[8],ia[3],mmid;
  double vol[8];
  int two[3]={2,2,2};
  
  name=readname(fpin);
  namevint=readname(fpin);
  oldnew=read1charname(fpin);
  printout("normal","volume integrate %s to %s, oldnew=%c\n",name,namevint,oldnew);
  
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
  v=(double *)need(name);
  
  vint=(double *)find(namevint);
  if (vint==0) 
  { 
    vint=(double *)createarray(namevint,noindfwpts[0],'d',0);
    Loop(i,0,noindfwpts[0]) vint[i]=0;
  }
  else if (oldnew=='n' || oldnew=='1')
    Loop(i,noindfwpts[its+1],noindfwpts[ite+2]) vint[i]=0;
  
  geom8init();
  geomcinit();
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* each cont. c.v. */
  {
    if (wherep[In4(ip,i4dp)]<0) continue; /* do only for valid volumes */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
    geom8vol(vol,ip);
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    {
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }
    Loop(k,0,8) vint[wherefw[mc[k]]]+=vol[k]*v[mmid];
  } /* each cont c.v. */
}	 
