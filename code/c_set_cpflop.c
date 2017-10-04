/* contains c_set_cpflop, coefprint1, omitpsleep */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
void coefprint1(int i, int *cn, int **ci, double **cc, int jrep, FILE *fprint); 
int omitpsleep(int iw, int *cn, int **ci, double **cc, int jrep, int *matchpc, double *cpsleep, int iallp); 

void c_set_cpflop(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*noindfwpts,*where,*matchpc; /* needed arrays */
  int *whois;
  double *cpda,*cpsleep;
  char *clt;
  int *itrange; double *roundoff; /* use if available */
  int *cpflop_n, **cpflop_i; double **cpflop_c; /* create or modify */
  
  int i,j,n,it,is,iw,ii[4],ia[4],ip[4],i4dp[4],*newai,nn,iwmin,iwmax,iallp,ipt[4];
  int its,ite,neqs,neqe;
  double tol=1.e-11;
  double *cnew;
  int itwo[3]={2,2,2};
  int nyexit=0;
  int nyprint=1;  /* set to 1 or 2 for debug prints */
  
  i4d=(int *)need("idim4d");      /* get or create needed arrays */
  Loop (i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  noindfwpts=(int *)need("noindfwpts");
  where=(int *)need("wherefw");
  whois=(int *)need("whoisfw");
  matchpc=(int *)need("matchpc");
  cpda=(double *)need("cpda"); 
  cpsleep=(double *)need("cpsleep"); 
  clt=(char *)need("clt");
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; neqs=0; neqe=noindfwpts[0]; }
  else {its=itrange[0]; ite=itrange[1]; neqs=noindfwpts[its+1]; neqe=noindfwpts[ite+2];}
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  cpflop_n=(int *)find("cpflop_n");
  n=noindfwpts[0];
  if (cpflop_n==0)
  {
    cpflop_n=(int *)createarray("cpflop_n",n,'i',0);
	 cpflop_i=(int **)createarray("cpflop_i",n,'p',0);
	 cpflop_c=(double **)createarray("cpflop_c",n,'p',0);
	 printout("normal","clear new cpflop\n");
	 Loop(i,0,n) {cpflop_n[i]=0; cpflop_i[i]=0; cpflop_c[i]=0;}
  }
  else
  { 
    cpflop_i=(int**)need("cpflop_i"); 
	 cpflop_c=(double**)need("cpflop_c"); 
    /* clear arrays to be reset, may need to consider nyrepeat when it is reformed*/
    /* printout("normal","clear for reset\n"); */
	 Loop(iw,neqs,neqe)
    if (cpflop_n[iw]>0) 
    {
      free(cpflop_i[iw]); free(cpflop_c[iw]);
      cpflop_n[iw]=0; cpflop_i[iw]=0;  cpflop_c[iw]=0;
    }
	 /*printout("normal","arrays cleared for reset\n"); */
  } /* not new  */
  /*printout("normal","ready to set cpflop\n"); */
  iwmin=noindfwpts[0]; iwmax=0;
  Loop(it,its,ite+1)   /* reset points */
  { 
    ii[3]=it; ip[3]=it;
    Loop3(ii,0,i4d)
    { 
      is=In4(ii,i4d); iw=where[is];
      if (iw<0)  continue;  /* do only for flow points */
      if (clt[is]=='w') continue;
      iwmin=min(iwmin,iw); iwmax=max(iwmax,iw);
      n=cpflop_n[iw];
      newai=(int *)smalloca(n+8,'i');
      cnew=(double *)smalloca((n+8)*3,'d');
      nn=n+8;
      if (n>0)
      {
        Loop(i,0,n) 
        {
          newai[i]=cpflop_i[iw][i]; 
          cnew[i]=cpflop_c[iw][i]; 
          cnew[i+nn]=cpflop_c[iw][i+n];
          cnew[i+2*nn]=cpflop_c[iw][i+2*n];
        }
        free(cpflop_i[iw]); free(cpflop_c[iw]); cpflop_i[iw]=0; cpflop_c[iw]=0;
      }
      Loop3(ia,0,itwo)
      {
        Loop (j,0,3) ip[j]=ii[j]+ia[j];
        newai[n]=In4(ip,i4dp);
        j=ia[0]+2*ia[1]+4*ia[2];
        cnew[n]=cpda[j+24*is];
        cnew[n+nn]=cpda[j+8+24*is];
        cnew[n+2*nn]=cpda[j+16+24*is];
        n++;
      }
      cpflop_n[iw]=nn;
      cpflop_i[iw]=newai;
      cpflop_c[iw]=cnew;
    } /* loop3 */
  } /* it */
  /* resolve sleeping points and consolidate  similar points */
  if (nyprint>0)  printout("normal","iw range %d %d\n",iwmin,iwmax);
  Loop(iw,iwmin,iwmax+1)   
  {
    if (nyprint>1) coefprint1(iw,cpflop_n,cpflop_i,cpflop_c,3,fprint);
    coefcombine(iw,cpflop_n,cpflop_i,cpflop_c,3,tol);  
    n=1;
    Loop(j,0,i4dp[0]+i4dp[1]+i4dp[2]) 
    {
      if (n>0)
      { 
        if (nyprint>1) coefprint1(iw,cpflop_n,cpflop_i,cpflop_c,3,fprint);
        n=omitpsleep(iw,cpflop_n,cpflop_i,cpflop_c,3,matchpc,cpsleep,iallp);
        if (n==0) break;
        if (nyprint>1) coefprint1(iw,cpflop_n,cpflop_i,cpflop_c,3,fprint);
        coefcombine(iw,cpflop_n,cpflop_i,cpflop_c,3,tol);  
      }
    }
    if (nyprint>1) coefprint1(iw,cpflop_n,cpflop_i,cpflop_c,3,fprint); 
    if (n>0) 
    {
      printout("error c_set_cpflop"," ERROR  unresolved sleeping points/n"); 
      coefprint1(iw,cpflop_n,cpflop_i,cpflop_c,3,fprint); 
      Loop(nn,0,cpflop_n[iw])
      {
        if (matchpc[cpflop_i[iw][nn]]>=0) 
          printout("error","pt %d matchpc %d %d\n",cpflop_i[iw][nn],
                  matchpc[cpflop_i[iw][nn]],matchpc[cpflop_i[iw][nn]+iallp]);
      }
      nyexit=1;
    }
  }
  if (nyexit==1) 
  {
    printout("error","exiting due to sleeping pts error\n"); 
    exitm4d(0); 
  }
  if (nyprint>0)
  { 
    double cmax,cmaxa=0,cr,crmax=0,sum[3];
    printout("normal","checking cpflop consistency\n");
    Loop(iw,iwmin,iwmax+1)
    {
      if (cpflop_n[iw]==0) continue;
      cmax=0;
      Loop(i,0,3) sum[i]=0;
      n=cpflop_n[iw];
      Loop(i,0,n) Loop(j,0,3)
      {
        cmax=max(cmax,abs(cpflop_c[iw][i+j*n]));
        sum[j]+=cpflop_c[iw][i+j*n];
      }
      if (cmax==0) continue;
      cmaxa=max(cmaxa,cmax);
      cr=0;
      Loop(j,0,3) cr=max(cr,abs(sum[j])/cmax);
      if (cr>.001 && cr*cmax/cmaxa> .00001)
      {
        iexpand(whois[iw],i4d,ipt);
        printout("warning c_set_cpflop","cpflop sum warning at iw=%d, ii %d %d %d %d, sum= %lg %lg %lg cmax= %lg\n",
                iw,ipt[0],ipt[1],ipt[2],ipt[3],sum[0],sum[1],sum[2],cmax);
      }
      crmax=max(crmax,cr);
    }
    printout("normal","    maximum sum/cmax = %lg\n",crmax);
    if (crmax>.001) 
    {
      printout("error c_set_cpflop","max sum/cmax > .001, check grid/bndry conditions, stop\n");
      exitm4d(0);
    }
  }
}
/* --------------------------- */
void coefprint1(int i, int *cn, int **ci, double **cc, int jrep, FILE *fprint)
{ 
  int j,n,k;
  n=cn[i];
  printout("normal"," point %d no of coef %d:  ic c1 c2 ...\n",i,n);
  Loop(k,0,n)
  {
    printout("normal","          %d ",ci[i][k]);
    Loop(j,0,jrep) printout("normal"," %lg",cc[i][k+j*n]);
    printout("normal","\n");
  }
}
/* --------------------------- */
int omitpsleep(int iw, int *cn, int **ci, double **cc, int jrep, int *matchpc, double *cpsleep, int iallp)
{ 
  int i,j,k,nn,np,iaddtot,nydo,*newai,*iadd,m1,m2; double *c,*cnew,*cmax;
  
  iaddtot=0; nydo=0;
  nn=cn[iw]; 
  if (nn>0)
  {
    iadd=(int*)smalloca(2*nn,'i');
    cmax=(double*)smalloca(nn,'d');
    c=cc[iw];
    Loop(j,0,nn)   /* first just count sleeping pts and no of added coef slots needed */
    { 
      i=ci[iw][j]; cmax[j]=0;
      if (matchpc[i]>0) /* then create cmax and check */
      { 
        Loop(k,0,3) cmax[j]=max(cmax[j],abs(c[j+nn*k]));
        if (cmax[j]>0)
        {
          nydo++;
          /*	   printout("normal","omit: eq %d pt %d matchpc %d %d\n",iw,i,matchpc[i],matchpc[i+iallp]); */
          Loop(k,0,nn) {if (matchpc[i]==ci[iw][k]) break;}
          if (k==nn) {iadd[iaddtot]=matchpc[i]; iaddtot++;}
          Loop(k,0,nn) {if (matchpc[i+iallp]==ci[iw][k]) break;}
          if (k==nn) {iadd[iaddtot]=matchpc[i+iallp]; iaddtot++;}
        }
      }
    }
    /*	  printout("normal","iaddtot %d \n",iaddtot); */
    if (iaddtot>0)  /* redo array if needed */
    {
      np=nn+iaddtot;
      newai=(int *)smalloca(np,'i');
      cnew=(double *)smalloca(np*3,'d');
      Loop(j,0,nn) 
      {
        newai[j]=ci[iw][j]; 
        Loop(k,0,3) cnew[j+np*k]=c[j+nn*k];
      }
      Loop(j,0,iaddtot)
      {
        newai[j+nn]=iadd[j];
        Loop(k,0,3) cnew[j+nn+np*k]=0;
      }
      free(ci[iw]); ci[iw]=newai;
      free(cc[iw]); cc[iw]=cnew;
      cn[iw]=np;
    }
    /*	   printout("normal","nydo %d\n",nydo); */
    if (nydo>0)  /* add in sleep coef */
    { 
      np=cn[iw]; 
      c=cc[iw];
      Loop(j,0,nn)  /* loop only over old points */
      { 
        i=ci[iw][j];
        if (cmax[j]>0)
        { 
          Loop(k,0,np) {if (matchpc[i]==ci[iw][k]) break;}
          m1=k; 
          Loop(k,0,np) {if (matchpc[i+iallp]==ci[iw][k]) break;}
          m2=k;
          if (m1==np || m2==np) printout("normal","error m1 m2 np %d %d %d\n",m1,2,np);
          Loop(k,0,3)
          { 
            c[m1+k*np]+=c[j+k*np]*cpsleep[i];
            c[m2+k*np]+=c[j+k*np]*cpsleep[i+iallp];
            c[j+k*np]=0;
          }
        }
      }
    }
    free(iadd); free(cmax);
  }
  return(nydo);
}
