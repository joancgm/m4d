/* contains c_momcam */
#include "global.h"
#include "eqnsubs.h"
/* abbreviated momentum coefficient cam for all 3 components, new formula*/

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_momcamddt(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherefw,*noindfwpts,*match,*wherep;  /* needed arrays */
  double *u[3],*rho;
  char *cltp;
  int *coef_n; double **coef_c; 
  int *itrange; double *zrotation; /* use if available */
  double *cama, *camb, *cplus;  /* set */
  int its,ite;  /* use itrange or all */
  
  int i,j,k,n,ip[4],i4dp[4],neqst,neqend,ia[3],ig[4],ipt,ieq;
  int two[3]={2,2,2};
  double cg[8],gradv[8][8][3],gradvu[8][3][3],rhom;
  double vol[8],fmid[8][3];
  double rot[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double cc,cmax;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  itrange=(int *)find("itrange");
  if (itrange>0) {its=itrange[0]; ite=itrange[1]; }
  else {its=0; ite=i4d[3]-1; }
  noindfwpts=(int *)need("noindfwpts");
  neqst=noindfwpts[its+1]; neqend=noindfwpts[ite+2];
  wherefw=(int *)need("wherefw");
  match=(int *)need("match");
  wherep=(int *)need("wherep");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  zrotation=(double *)find("zrotation");
  if (zrotation>0) 
  {  
    rot[0][1]=-zrotation[0]; 
    rot[1][0]=zrotation[0];
  }
  coef_n=(int *)need("coef_n");
  coef_c=(double **)need("coef_c");
  cltp=(char *)need("cltp");
  geom8init();  
  geomcinit();
  
  cplus=(double *)createarray("cplus",noindfwpts[0],'d',0);
  eqncplus(neqst,neqend,coef_n,coef_c, cplus);   /* + coefficients */
  
  cama=(double *)createarray("cama",noindfwpts[0],'d',0);
  camb=(double *)createarray("camb",noindfwpts[0],'d',0);
  Loop(i,neqst,neqend) {cama[i]=0; camb[i]=0;} /* clear before reseting */
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* loop over all possible cont c.v. */
  {
    if (wherep[In4(ip,i4dp)]<0) continue; 
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    ig[3]=ip[3];
    geom8gradvol(gradv,ip); 
    geom8volfmid(vol,fmid,ip);
    geomcvave(cg,ip); 
    rhom=0;
    Loop(i,0,8) Loop(j,0,3) Loop(k,0,3) gradvu[i][j][k]=0;
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
      ipt=In4(ig,i4d);
      j=ia[0]+2*ia[1]+4*ia[2];
      rhom+=cg[j]*rho[ipt];
      Loop(i,0,8) Loop(k,0,3) Loop(n,0,3) gradvu[i][k][n]+=gradv[i][j][k]*u[n][ipt];
    }
    /* Loop(i,0,8) fixdudx(gradvu[i]); */
    
    Loop3(ia,0,two)
    {
      Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
      ieq=wherefw[In4(ig,i4d)];
      j=ia[0]+2*ia[1]+4*ia[2];
      if (vol[j]<=0) continue;
      if (ieq<0) continue;
      /* totsafe based on eqs 23,24 25 */
      cmax=0;
      Loop(i,0,3) 
      { cc=0;
        Loop(k,0,3) cc+=(gradvu[j][i][k]+2*vol[j]*rot[i][k])*(gradvu[j][k][i]+2*vol[j]*rot[k][i]);
        if (cc>0)  cc=sqrt(abs(cc));
        else cc=0; 
        cmax=max(cc,cmax);
      }
      cama[ieq]+=rhom*cmax;
      Loop(i,0,3) cmax=max(cmax,abs(gradvu[j][i][i]));  /* eq 19 */
      camb[ieq]+=rhom*cmax;
    }
  }
}
