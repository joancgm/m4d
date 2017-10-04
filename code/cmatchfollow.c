/* contains cmatchfollow */
#include "global.h"
/*     follow cmatch to valid point  */
#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))


int 	cmatchfollow(int kin, char *cmatch, int *wherepx, int ippt, int *i4dp, 
                   int iallp, int nyuser, int *iptout, int *kout)
{  
  int i,k,kop,knew,ipt,iptnew,iter,itermax,ip[4],L,side,di,nythrurepeat=0;
  char cm, cx[9]="iIjJkKrR";
  int bug=0;
  
  *iptout=ippt;
  *kout=kin;
  ipt=ippt;
  k=kin;
  iexpand(ipt,i4dp,ip);
  itermax=i4dp[0]+i4dp[1]+i4dp[2];
  
  Loop(iter,0,itermax)
  { 
    if (wherepx[ipt]>=0) break;
    L=k/2; side=k%2; di=2*side-1;
    ip[L]+=di; 
    iptnew=In4(ip,i4dp);
    kop=2*L+1-side;
    cm=cmatch[iptnew+kop*iallp];
    knew=-1; 
    Loop(i,0,6) {if (cm==cx[i]) {knew=i; break;}}
    if (bug>0 && iter==0) 
    {
      printout("normal","cmatchfollow ipt %d ip %d %d %d %d k %d cm %c L %d side %d di %d kop %d cmp %c knew %d\n",
             ipt,ip[0],ip[1],ip[2],ip[3],k,cmatch[ipt+k*iallp],L,side,di,kop,cm,knew);
    }
    
    if (ip[L]==0 || ip[L]==i4dp[L]-1) /* check if off edge */
    { 
      if (nyuser==1 && (cm=='r' || cm=='R'))  /* check repeat follow */
      { 
        if (ip[L]==0) {ip[L]=i4dp[L]-2; nythrurepeat=-1; }
        else { ip[L]=1; nythrurepeat=1; }
        ipt=In4(ip,i4dp);
        continue;
      }
      else
      {
        ipt=-ipt; break;
      }
    }
    ipt=iptnew;
    if (wherepx[ipt]>=0) break;  /* check wherepx */
    if (iter==itermax-1) {ipt=-ipt; break;} /* find matching side ready for next point */
    if (knew==-1) {ipt=-ipt; break;}
    k=knew;
  }
  
  if (bug>0) printout("normal","        iptout %d  ip %d %d %d %d kout %d\n",ipt,ip[0],ip[1],ip[2],ip[3],k);
  *iptout=ipt; *kout=k;
  return (nythrurepeat);		
}
