/* contains c_copy_mtop */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_copy_mtop(FILE *fpin, FILE *fprint)   
/*  copy mid-point to p-point for arrays of type double */
{ 
  char *namefrom, *nameto; int nysym;  /* if nysym=0 outside pts set to zero */
  int *i4d; char *csym;
  int i4dp[4],i4dm[4],i,j,k,n,m,ip[4],im[4],nyout,iallp,iallm;
  double *vf,*vto;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iallp=Prod4(i4dp);
  
  namefrom=readname(fpin);
  nameto=readname(fpin);
  nysym=readint(fpin);
  printout("normal","copy %s to %s nysym= %d\n",namefrom,nameto,nysym);
  i=arraysize(namefrom);
  if (i != iallm) 
  {
    printout("error c_copy_mtop"," array %s has size %d not %d, cannot copy\n",namefrom,i,iallm); 
    exitm4d(0); 
  }
  vf=(double *)need(namefrom);
  vto=(double *)createarray(nameto,iallp,'d',0);
  
  Loop(ip[3],0,i4dp[3]) Loop3(ip,0,i4dp)
  { 
    Loop(i,0,3) im[i]=max(0,min(i4dm[i]-1,ip[i]-1)); 
    im[3]=ip[3];
    nyout=0; Loop(i,0,3) if (im[i]!=ip[i]-1) nyout=1;
    if (nyout == 1 && nysym == 0) vto[In4(ip,i4dp)]=0;
    else vto[In4(ip,i4dp)]=vf[In4(im,i4dm)];
  }
  if (nysym>0)
  { 
    csym=(char *)need("csym");
    Loop(i,0,3) 
    {  
      j=(i+1)%3; k=(i+2)%3;
      if (csym[6+2*i] == 'r')
        Loop(ip[3],0,i4dp[3]) Loop(ip[j],0,i4dp[j]) Loop(ip[k],0,i4dp[k])
      { 
        ip[i]=0; m=In4(ip,i4dp);
        ip[i]=i4dp[i]-1; n=In4(ip,i4dp); 
        vto[n]=.5*(vto[m]+vto[n]);
        vto[m]=vto[n];
      }
    }
  }
}
