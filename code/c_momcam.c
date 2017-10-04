/* contains c_momcam */
#include "global.h"
#include "eqnsubs.h"
/* abbreviated momentum coefficient cam for all 3 components */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_momcam(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherefw,*noindfwpts,*match,*wherep,*itrange;  /* needed arrays */
  double *u[3],*rho;
  char *cltp;
  int *coef_n; double **coef_c; 
  double *cam, *cplus;  /* set */
  double relax;
  
  int its,ite;  /* use itrange or all */
  
  int i,j,k,n,ip[4],i4dp[4],neqst,neqend,ia[3],ig[4],ipt,ieq;
  int two[3]={2,2,2};
  double cg[8],gradv[8][8][3],gradvu[8][3][3],rhom,*camold=0,camt[3];
  
  relax=readdouble(fpin);  /* relaxation with old value */
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
  cltp=(char *)need("cltp");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  geom8init();  
  geomcinit();
  coef_n=(int *)need("coef_n");
  coef_c=(double **)need("coef_c");
  
  if (arraysize("cam")!=noindfwpts[0]) relax=0;
  printout("normal","cami=max of caminew and  %lg caminew +%lg camiold\n",1.-relax,relax);
  
  cplus=(double *)createarray("cplus",noindfwpts[0],'d',0);
  eqncplus(neqst,neqend,coef_n,coef_c, cplus);   /* + coefficients */
  cam=(double *)createarray("cam",noindfwpts[0],'d',0);
  
  if (relax>0)
  { 
    camold=(double*)tmalloca(noindfwpts[0],'d');
    Loop(i,neqst,neqend) camold[i]=cam[i];
  }
  Loop(i,neqst,neqend) cam[i]=0; /* clear before ressing */
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* loop over all possible cont c.v. */
  {
    if (wherep[In4(ip,i4dp)]<0) continue; 
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    ig[3]=ip[3];
    geom8gradvol(gradv,ip); 
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
    Loop(i,0,8) Loop(j,0,3) Loop(k,0,3) gradvu[i][j][k]=sqrt(abs(gradvu[i][j][k]));
    
    Loop3(ia,0,two)
    {
      Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
      ieq=wherefw[In4(ig,i4d)];
      j=ia[0]+2*ia[1]+4*ia[2];
      if (ieq<0) continue;
      Loop(k,0,3) camt[k]=rhom*(gradvu[j][k][0]*gradvu[j][0][k]
                                +gradvu[j][k][1]*gradvu[j][1][k]
                                +gradvu[j][k][2]*gradvu[j][2][k]);
      camt[0]=max(camt[0],max(camt[1],camt[2]));
      cam[ieq]+=camt[0];				
    }
  }
  if (relax>0)
    Loop(i,neqst,neqend) 
  { 
    camold[i]=relax*camold[i]+(1.-relax)*cam[i];
    if (camold[i]>cam[i]) cam[i]=camold[i];
  }
}


