/* contains c_coefconv, c_coefconvstep, coefconv */
#include "global.h"
/* calculate coefficient for convection term  based on rhoU grad prop dVol */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
void coefconv(int nylin,int jconv, FILE *fprint);

/* ----------------------------------- */
void c_coefconv(FILE *fpin, FILE *fprint)
{ 
  int nylin=1;
  int jconv;   /* input parameter */
  jconv=readint(fpin);
  printout("normal"," convection=coef %d",jconv);
  coefconv(nylin,jconv,fprint);
}
/* ----------------------------------- */
void c_coefconvstep(FILE *fpin, FILE *fprint)
{ 
  int nylin=0;
  int jconv;   /* input parameter */
  jconv=readint(fpin);
  printout("normal"," convection=coef %d",jconv);
  coefconv(nylin,jconv,fprint);
}
/* ----------------------------------- */
void coefconv(int nylin, int jconv, FILE *fprint)
{ 
  int *i4d,*wherep,*noindfwpts,*wherefw,*match,*whoisfw,*nocoefs;     /* needed arrays */
  double *u[3],*rho;
  char *cltp;
  int *itrange; double *roundoff; /* use if available */
  int *coef_n, **coef_i; double **coef_c;  /* modify */
  int its,ite;
  int i,j,k,ia[3],ib[3],ic[4],ig[4],ip[4],i4dp[4],ieq,ipt,iptc[8],iw[8];
  int two[3]={2,2,2};
  double cg[8],rhou[3],gradv[8][8][3],coefa[8][8];
  double ddI[3][3][2];  /* for stepwide discretization */
  double tol=1.e-8;
  
  nocoefs=(int *)need("nocoefs");
  if (jconv>nocoefs[0]-1)
  { 
    printout("error coefconv"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
    exitm4d(0);
  }
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  printout("normal",", set for its=%d to ite=%d\n",its,ite);
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  noindfwpts=(int *)need("noindfwpts");
  wherefw=(int *)need("wherefw");
  whoisfw=(int *)need("whoisfw");
  match=(int *)need("match");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  geom8init();  
  geomcinit();
  coef_n=(int *)need("coef_n");
  coef_i=(int **)need("coef_i");
  coef_c=(double **)need("coef_c");
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])  /* clear parts of array to be redone */
  if (coef_n[i]>0) 
  { Loop(j,0,coef_n[i]) coef_c[i][j+coef_n[i]*jconv]=0; }
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* loop over all possible cont c.v. */
  {    
    ic[3]=ip[3];  ig[3]=ip[3];   
    /* do if valid cont c.v. */
    if (wherep[In4(ip,i4dp)]>=0 && cltp[In4(ip,i4dp)]!='s' && cltp[In4(ip,i4dp)]!='S')     
    { 
      geomcvave(cg,ip);  /*   set rhou    */
      /*	  printout("normal","ip %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]); */
      Loop(i,0,3) rhou[i]=0;
      Loop(i,0,3) ic[i]=ip[i]-1;
      Loop3(ia,0,two)
      { 
        Loop(i,0,3) ig[i]=ic[i]+ia[i];
        ipt=In4(ig,i4d);
        j=ia[0]+2*ia[1]+4*ia[2];	
        Loop(i,0,3) rhou[i]+=cg[j]*rho[ipt]*u[i][ipt];
      }
      /*	  printout("normal","rhou %lg %lg %lg\n",rhou[0],rhou[1],rhou[2]); */
      if (Sqr(rhou)>0)  /* non-zero velocity */
      {   
        geom8gradvol(gradv,ip);  /* 1/8th gradients */
        if (nylin==0) /* fold in for stepwise */
        {  
          Loop3(ia,0,two) /* eight corner volumes */
          { 
            j=ia[0]+2*ia[1]+4*ia[2];
            Loop(i,0,3) Loop(k,0,3) { ddI[i][k][0]=0; ddI[i][k][1]=0; }
            Loop3(ib,0,two) /* eight corner points */
            {
              Loop(i,0,3) Loop(k,0,3) 
              ddI[i][k][ib[i]]+=gradv[j][ib[0]+2*ib[1]+4*ib[2]][k];
            }
            Loop(i,0,8) Loop(k,0,3) gradv[j][i][k]=0;
            Loop(k,0,3) Loop(i,0,2)
            {  
              gradv[j][i+2*ia[1]+4*ia[2]][k]+=ddI[0][k][i];
              gradv[j][ia[0]+2*i+4*ia[2]][k]+=ddI[1][k][i];
              gradv[j][ia[0]+2*ia[1]+4*i][k]+=ddI[2][k][i];
            }
          }
        }  /* end stepwise fold */
        
        Loop(i,0,8) Loop(j,0,8) coefa[i][j]=0;
        Loop3(ia,0,two)  /* eight corners */
        { 
          j=ia[0]+2*ia[1]+4*ia[2];
          Loop(i,0,3) ig[i]=ic[i]+ia[i];
          iptc[j]=match[In4(ig,i4d)];
          /*	  printout("normal","ia %d %d %d ig %d %d %d %d\n",ia[0],ia[1],ia[2],ig[0],ig[1],ig[2],ig[3]); */
          iw[j]=wherefw[iptc[j]];
          if (iw[j]<0) continue;
          Loop3(ib,0,two)
          { 
            Loop(i,0,3) ig[i]=ic[i]+ib[i];
            k=ib[0]+2*ib[1]+4*ib[2];
            coefa[j][k]=Dot(rhou,gradv[j][k]);
          }
        } /* ia eight corners */
        Loop(i,0,8) coefcadd(iw[i],coef_n,coef_i,coef_c,nocoefs[0],jconv,8,iptc,coefa[i]);
      } /* non-zero velocity */
    } /* valid cont c.v. */
  } /* ip loop */  
  Loop(ieq,noindfwpts[its+1],noindfwpts[ite+2]) 
  {
    coefcombine(ieq,coef_n,coef_i,coef_c,nocoefs[0],tol);
    coeforder(whoisfw[ieq],ieq,coef_n,coef_i,coef_c,nocoefs[0]);  
  }
}

