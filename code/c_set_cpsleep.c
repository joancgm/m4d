/* contains c_set_cpsleep */
#include "global.h"

/* coefficients for sleeping points */
/* set  cpsleep */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

void c_set_cpsleep(FILE *fpin, FILE *fprint)
{ 
  int *i4d, *matchpc; char *csym;  double *xpc; /* needed arrays */
  double  *cpsleep; /* arrays created */
  
  int i,n,k,k1,k2,ka,kb,i4dp[4],ip[4],iallp,L,ipma[4],ipmb[4],ipmc[4];
  double xa[3],xb[3],dxr[3],f;
  
  i4d=(int *)need("idim4d");
  csym=(char *)need("csym");
  matchpc=(int *)need("matchpc");
  xpc=(double *)need("xyzp");
  
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  
  cpsleep=(double *)createarray("cpsleep",iallp*2,'d',0);
  Loop(i,0,iallp*2) cpsleep[i]=0;
  
  Loop(ip[3],0,i4dp[3])  Loop3(ip,0,i4dp)
  {
    k=In4(ip,i4dp);
    if (matchpc[k] <0) continue;
    ka=matchpc[k];
    kb=matchpc[k+iallp];
    /*printout("normal","k %d ka %d kb %d\n",k,ka,kb); */
    cpsleep[k]=1;
    cpsleep[k+iallp]=0;
    if (ka!=kb) /* set spatial interp coef  based on projected distance*/
    { 
      Loop(i,0,3)
      { 
        xa[i]=xpc[ka+iallp*i]-xpc[k+iallp*i];
        xb[i]=xpc[kb+iallp*i]-xpc[k+iallp*i];
      }
      iexpand(ka,i4dp,ipma);
      iexpand(kb,i4dp,ipmb);
      n=0;
      Loop (L,0,3)     /* check repeat bndry dx */
      if (csym[6+2*L]=='r' && (ipma[L]-ip[L])*(ipmb[L]-ip[L])>0)  
      { 
        n=1;
        Loop(i,0,4) ipmc[i]=ip[i];
        ipmc[L]=0; k1=In4(ipmc,i4dp);
        ipmc[L]=i4dp[L]-1; k2=In4(ipmc,i4dp);
        /*printout("normal","L %d k %d k1 %d k2 %d ka %d kb %d \n",L,k,k1,k2,ka,kb); */
        Loop(i,0,3) dxr[i]=xpc[k2+iallp*i]-xpc[k1+iallp*i];
        /*printout("normal","xa  %lg %lg %lg\n",xa[0],xa[1],xa[2]);
         printout("normal","xb  %lg %lg %lg\n",xb[0],xb[1],xb[2]);
         printout("normal","dxr %lg %lg %lg\n",dxr[0],dxr[1],dxr[2]); */
        if (abs(ipma[L]-ip[L])>abs(ipmb[L]-ip[L]))
        { 
          if (ipma[L]-ip[L]>0) Loop(i,0,3) xa[i]-=dxr[i];
          else Loop(i,0,3) xa[i]+=dxr[i];
          /*printout("normal","xa  %lg %lg %lg\n",xa[0],xa[1],xa[2]); */
        }
        else 
        { if (ipmb[L]-ip[L]>0) Loop(i,0,3) xb[i]-=dxr[i];
        else Loop(i,0,3) xb[i]+=dxr[i];
          /*printout("normal","xb  %lg %lg %lg\n",xb[0],xb[1],xb[2]); */
        }
      }
      Loop(i,0,3) xb[i]=xa[i]-xb[i];
      f=Dot(xb,xb);  
      if (f>0) f=Dot(xa,xb)/f;
      /*if (n==1) printout("normal","f %lg\n",f); */
      f=max(0,min(1,f));
      cpsleep[k]=1-f;
      cpsleep[k+iallp]=f;
    }
  }
}



