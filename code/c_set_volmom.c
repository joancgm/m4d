/* contains c_set_volmom */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_set_volmom(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*noindfwpts,*wherefw,*matchpc; /* neeeded arrays */
  char *cltp;
  double *volmom; /* set  volmom */
  
  int i,nn,ip[4],i4dp[4],ia[3],ifpt[4],npt,ipp;
  int two[3]={2,2,2};
  double vol8[8];
  
  i4d=(int *)need("idim4d");      /* get needed arrays */
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  noindfwpts=(int *)need("noindfwpts");
  matchpc=(int *)need("matchpc");
  wherefw=(int *)need("wherefw");
  cltp=(char *)need("cltp");
  nn=noindfwpts[0];
  geom8init();   /* get needed arrays for vol8th */
  
  volmom=(double *)createarray("volmom",nn,'d',0);  /* create volmom */
  Loop(i,0,nn) volmom[i]=0;
  Loop(ip[3],0,i4d[3]) Loop3(ip,1,i4d)  /* loop over internal  p points */
  {
    ipp=In4(ip,i4dp);
    if (matchpc[ipp] >=0) continue;  /* do only for non-zero volumes */
    if (cltp[ipp]=='s' || cltp[ipp]=='S') continue;  /* omit solids */
    geom8vol(vol8,ip);
    Loop3(ia,0,two)   /* each "1/8th" volume */
    {
      Loop(i,0,3) ifpt[i]=ip[i]-1+ia[i]; ifpt[3]=ip[3];  /* flow pt index */
      npt=wherefw[In4(ifpt,i4d)];
      if (npt>=0)  /* point may be used */
        volmom[npt]+=vol8[ia[0]+2*ia[1]+4*ia[2]];
    } /* Loop  ia */
  }  /* loop ip */
}
