/* contains psleepcalc */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* set property p for sleeping p-points using cpsleep interpolation
 call with cpsleep for actual p points 
 and cpsleepm for center located i.e. mid point arrays */

void psleepcalc(FILE *fprint, double *p, double *cpsleep)
{ 
  int *i4d,*noppts,*noindppts,*matchpc;  /* needed arrays */
  char *cltp;
  int *itrange; /* use if available */
  int i,j,k,i4dp[4],iallp;
  double c[22],sum; int ipa[22],n,ny;
  int is,ie,iperr[4],nn,kk,ifac[3];
  int bug=0; double csum;
  
  i4d=(int *)need("idim4d");      /* get or create needed arrays */
  Loop (i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  ifac[0]=1; ifac[1]=i4dp[0]; ifac[2]=i4dp[0]*i4dp[1];
  iallp=Prod4(i4dp);
  noppts=(int *)need("noppts");
  noindppts=(int *)need("noindppts");
  matchpc=(int *)need("matchpc");
  cltp=(char *)need("cltp");
  is=0; ie=iallp;
  itrange=(int *)find("itrange");
  if (itrange>0) 
  {
    is=itrange[0]*i4dp[0]*i4dp[1]*i4dp[2];
    ie=(itrange[0]+1)*i4dp[0]*i4dp[1]*i4dp[2];
  }
  
  Loop(i,is,ie)
  if (matchpc[i]>=0)
  {
    n=1; c[0]=1.; ipa[0]=i; ny=1;
    while (ny)
    {
      ny=0;
      Loop(j,0,n)
      if (c[j]>0)
        if (matchpc[ipa[j]]>=0)
        {
          ny=1;
          c[n]=c[j]*cpsleep[ipa[j]];
          ipa[n]=matchpc[ipa[j]];
          c[n+1]=c[j]*cpsleep[ipa[j]+iallp];
          ipa[n+1]=matchpc[ipa[j]+iallp];
          n=n+2;
          c[j]=0;
          if (n+2>20) 
          { 
            iexpand(i,i4dp,iperr);
            printout("warning psleepcalc"," psleep error %d  %d %d %d %d matchpc %d %d\n",
                    i, iperr[0],iperr[1],iperr[2],iperr[3],matchpc[i],matchpc[i+iallp]); 
            Loop(nn,0,n)
            { 
              kk=ipa[nn]; iexpand(kk,i4dp,iperr);
              printout("warning","nn %d ,using pt %d %d %d %d %d  c %lg\n",
                      nn,kk, iperr[0],iperr[1],iperr[2],iperr[3],c[nn]); 	   
            }
            ny=0; 
            break;
          }
        }
    }
    sum=0;
    Loop(j,0,n) sum+=c[j]*p[ipa[j]];
    p[i]=sum;
    if (bug>1)
    { 
      printout("normal"," psleep at %d = %lg",i,sum); 
      iexpand(i,i4dp,iperr);
      printout("normal","   %d %d %d %d mat %d %d cps %lg %lg\n",
              iperr[0],iperr[1],iperr[2],iperr[3],matchpc[i],matchpc[i+iallp],cpsleep[i],cpsleep[i+iallp]); 
      sum=0; csum=0;
      Loop(j,0,n)
      { 
        sum+=c[j]*p[ipa[j]];
        csum+=c[j];
        iexpand(ipa[j],i4dp,iperr);
        printout("normal"," using p at %d   %d %d %d %d =%lg times %lg  sum %lg csum %lg mat %d %d cps %lg %lg\n",
                ipa[j], iperr[0],iperr[1],iperr[2],iperr[3],p[ipa[j]],c[j],sum,csum,
                matchpc[ipa[j]],matchpc[ipa[j]+iallp],cpsleep[ipa[j]],cpsleep[ipa[j]+iallp]); 
      }
    }
  }
  Loop(i,is,ie) if (cltp[i]=='S') /* set inside solid based on near wall values */
  { 
    iexpand(i,i4dp,iperr);
    if (bug>1) printout("normal"," p internal at %d   %d %d %d %d",i,iperr[0],iperr[1],iperr[2],iperr[3]); 
    n=0;
    Loop(j,0,3) for (k=-1;k<2;k+=2)
    { 
      kk=i+k*ifac[j];
      if (kk<0 || kk>=iallp) continue;
      if (cltp[kk]!='F') continue;
      if (n==0) p[i]=0;
      p[i]+=p[kk]; n++;
      if (bug>1) 
      {
        iexpand(kk,i4dp,iperr);
        printout("normal"," using p at %d   %d %d %d %d =%lg\n",kk,iperr[0],iperr[1],iperr[2],iperr[3],p[kk]);
      }
    }
    if (n>0) p[i]/=(double)n;
    if (bug>1) printout("normal"," p=%lg\n",p[i]);							  
  }
}
