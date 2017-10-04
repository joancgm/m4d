/* contains psolve1d */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* 1dsolver for pressure update dp,  i j or k direction */
/* update dp */
void psolve1d(int m, int *cn, int **ci, double **cc, int *i4dp, 
              int *wherep, char *cltp, int *ipfix, double *rhsc, double *dp)
{  
  double *dx,*cpn,*bbn,*ccn,*cmn,d;  /* for tdma */
  int ip[4],i,j,ifix=-1,ipta,ieq,m2,m3,ipt,ipp[4],max,ipfixi=-1;
  
  max=i4dp[m];
  dx=(double*)smalloc(max*5*sizeof(double));
  Loop(i,0,max*5) dx[i]=0;
  cpn=dx+max; bbn=cpn+max; ccn=bbn+max; cmn=ccn+max;
  
  ip[3]=0;
  m2=(m+1)%3; m3=(m+2)%3;
  /* coef and rhs */	
  Loop(ip[m],0,max) 
  {	
    Loop(ip[m2],0,i4dp[m2]) Loop(ip[m3], 0,i4dp[m3])
    { 
      if (ip[m]==ifix || ip[m]==ipfixi) continue;
      ipta=In4(ip,i4dp);
      ieq=wherep[ipta];
      if (ieq<0) continue;
      if (cn[ieq]==0) continue;
      bbn[ip[m]]+=rhsc[ieq];
      if (cltp[ipta]=='o')
      { 
        ifix=ip[m]; 
        ccn[ip[m]]=1; 
        cmn[ip[m]]=0; 
        cpn[ip[m]]=0; 
        bbn[ip[m]]=0;
      }
      if (ipfix>0 && ipfix[m]==ip[m])
      {
        ipfixi=ip[m]; 
        ccn[ip[m]]=1; 
        cmn[ip[m]]=0; 
        cpn[ip[m]]=0; 
        bbn[ip[m]]=0;
      }
      Loop(i,0,cn[ieq])
      { 
        ipt=ci[ieq][i];
        iexpand(ipt,i4dp,ipp);
        if (ipp[m]==ip[m]) ccn[ip[m]]+=cc[ieq][i];
        else if (ipp[m]==ip[m]-1) cmn[ip[m]]+=cc[ieq][i];
        else if (ipp[m]==ip[m]+1) cpn[ip[m]]+=cc[ieq][i];
        else if (cc[ieq][i]>0)  ccn[ip[m]]+=cc[ieq][i];
      }
    }
    if (ccn[ip[m]]==0) ccn[ip[m]]=1;
  }
  /* solve */				
  
  if (max==1) dx[0]=bbn[0]/ccn[0];
  else if (max==2)
  {
    dx[1] = (bbn[1]-cmn[1]*bbn[0]/ccn[0])/(ccn[1]-cmn[1]*cpn[0]/ccn[0]);
    dx[0] = bbn[0]/ccn[0] -dx[1]*cpn[0]/ccn[0];
  }
  else   /* max>2 */
  {
    cpn[0] = -cpn[0]/ccn[0];
    bbn[0] = bbn[0]/ccn[0];
    for (j=1;j<max-1;j++)
    { 
      d = ccn[j]+cmn[j]*cpn[j-1];
      bbn[j] = (bbn[j]-cmn[j]*bbn[j-1])/d;
      cpn[j] = -cpn[j]/d;
    }
    dx[max-1] = (bbn[max-1]-cmn[max-1]*bbn[max-2])
    /(ccn[max-1]+cmn[max-1]*cpn[max-2]);
    for (j=max-2;j>=0;j--) dx[j] = bbn[j]+cpn[j]*dx[j+1];
  }
  /* update */
  Loop(ip[m],0,max) 
  Loop(ip[m2],0,i4dp[m2]) Loop(ip[m3],0,i4dp[m3])
  {
    ipta=In4(ip,i4dp); 
    ieq=wherep[ipta];
    if (ieq<0) continue;
    if (cn[ieq]==0) continue;
    dp[ipta]+=dx[ip[m]];
  }
  free(dx);
}
