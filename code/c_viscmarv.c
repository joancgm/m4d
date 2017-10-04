/* contains c_viscmarv */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* calculate MARV model turbulent viscosity for q and bij 
 vmarv(6 ind. components)=.0735 rho sqrt(uiuj) q/om
 (for om equation multiply by .1/.0735= 1.36054 
 or for chanpipe consistency: 1.36666)
 */

void c_viscmarv(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep;   /* need */
  double *bij, *q, *om, *rho; 
  int *itrange; /* use if available */ 
  double *vmarv;  /* set */
  double prinstress(double (*re)[3], double (*retr)[3], double (*tr)[3], double tol);
  
  int i,k,i4dm[4],iallm,ip[4],i4dp[4],its,ite,ii[4],ia[3],mc[8],mmid,iall;
  double bm[6],qm,omm,rhom,cave[8],fac,vijdq[3],trd3;
  double uiujdk[3][3],uiujpdk[3][3],tr[3][3];
  
  int iprt=10;
  int ij2n[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
  int n2ii[6]    = {0,1,2,0,0,1};
  int n2jj[6]    = {0,1,2,1,2,2};
  int two[3]={2,2,2};
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallm=Prod4(i4dm);
  
  rho=(double *)need("rho");
  q=(double *)need("qturb");
  om=(double *)need("omturb");
  bij=(double *)need("bij");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1];}
  else {its=0; ite=i4d[3]-1; }
  geomcinit();
  vmarv=(double *) find("vmarv");
  if (vmarv==0)
  { 
    vmarv=(double*)createarray("vmarv",iallm*6,'d',0);
    Loop(i,0,6) vmarv[i]=0;
  }
  
  /* loop over contiuity control volumes */
  Loop(ip[3],its,ite+1) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));  /* viscosity index */
    Loop(i,0,6) vmarv[mmid+iallm*i]=0;   /* initialize as zero */
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }  
    if (wherep[In4(ip,i4dp)]<0) Loop(i,0,8) cave[i]=.125; /* simple ave */
    else geomcvave(cave,ip);   /* real cont. c.v. ave */
    Loop(i,0,6) bm[i]=0; qm=0; omm=0; rhom=0;
    Loop(k,0,8)
    {
      rhom+=cave[k]*rho[mc[k]];
      /* omm+=cave[k]/om[mc[k]]; */ /* assume 1/om varies linearly */
      omm+=cave[k]*log(om[mc[k]]);
      qm+=cave[k]*q[mc[k]];
      Loop(i,0,6)  bm[i] += cave[k]*bij[mc[k]+iall*i];
    }
    /* omm=1./omm; */
    omm=exp(omm);  
    if (qm<=0) continue;  /* can't set, leave as zero */
    Loop(i,0,3) Loop (k,0,3) uiujdk[i][k]=2*bm[ij2n[i][k]];
    Loop(i,0,3) uiujdk[i][i]+=2./3.;
    fac=prinstress(uiujdk,uiujpdk,tr,1.e-6);
    Loop(i,0,3) vijdq[i]=uiujpdk[i][i];
    if (vijdq[0]<0 || vijdq[1]<0 || vijdq[2]<0)  /* fix for relizability before sqrt */
    {
      i=0; if (vijdq[1]<vijdq[0]) i=1; if (vijdq[2]<vijdq[i]) i=2;
      trd3=(vijdq[0]+vijdq[1]+vijdq[2])/3.;
      fac=-trd3/(vijdq[i]-trd3);
      if (iprt>0)
        printout("warning c_viscmarv","non-realizable bij at ip: %d %d %d %d, fac=%lg used\n",
                ip[0],ip[1],ip[2],ip[3],fac);
      iprt--;
      Loop(i,0,3) vijdq[i]=(vijdq[i]-trd3)*fac+trd3;
    }
    Loop(i,0,3) vijdq[i]=sqrt(vijdq[i]);
    
    fac=0.0735*rhom*qm*qm/omm;
    Loop(k,0,6) vmarv[mmid+k*iallm]=fac*(tr[0][n2ii[k]]*tr[0][n2jj[k]]*vijdq[0]
                                         +tr[1][n2ii[k]]*tr[1][n2jj[k]]*vijdq[1]
                                         +tr[2][n2ii[k]]*tr[2][n2jj[k]]*vijdq[2]);
  }
  if (iprt<10) printout("warning c_viscmarv","realizability correction applied at %d points\n",-iprt-10);
}
