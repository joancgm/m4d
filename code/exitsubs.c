/* contains c_exitinit c_contcpcexit, c_contrhsexit */
#include "global.h"
/* subs for flow out pressure boundary conditions */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
void coefprint1(int i, int *cn, int **ci, double **cc, int jrep, FILE *fprint); 

/* set up for multiple  outflow regions with 
 flowout[0]= nout, number of regions
 then for each region: type  dL  ist[4] iend[4] ifix[4, flowoutprop.
 flowoutprop[nout]: type=1,2 pressure, type 3,4 flow rate
 type 1 fixed p at ifix (dp=0), flowoutprop ignored, dp/diL=uniform
 ( type 2 fixed pave, dp/diL=uniform , not yet coded)
 type 3 fixed flow rate,  dp/diL=uniform
 type 4 fixed flow rate, dp uniform 
 */
/* fixed pressure at each point does not need coding, it is the default */

/* ------------------------- */
void c_exitinit(FILE *fpin, FILE *fprint)
{
  int *i4d,*matchpc; /* need */
  int *flowout; double *flowoutprop;  /* create */
  int nout,*type,*dL,*ist,*iend,*ifix;
  
  int i,j,L,n,i4dp[4],ii[4], ipt;
  
  i4d=(int *)need("idim4d"); 
  matchpc=(int *)need("matchpc");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  nout=readint(fpin); 
  printout("normal"," %d outflow regions regions\n",nout);
  flowout=(int *)createarray("flowout",1+16*nout,'i',0);   
  flowoutprop=(double *)createarray("flowoutprop",nout,'d',0);   
  
  flowout[0]=nout;    /* read number of regions */
  n=1;
  Loop(i,0,nout)       /* read info for each region */
  {
    type=flowout+n; dL=type+1; ist=dL+1; iend=ist+4; ifix=iend+4; 
    type[0]=readint(fpin);
    dL[0]=readint(fpin);
    printout("normal","region %d, type %d, dL %d, ",i+1,type[0],dL[0]);
    Loop(j,0,4) 
    {
      ist[j]=readint(fpin)-1; 
      if (ist[j]<0) ist[j]=0; 
      if (ist[j]>i4dp[j]-1) ist[j]=i4dp[j]-1; 
    }
    printout("normal"," is[] %d %d %d %d, ",ist[0],ist[1],ist[2],ist[3]);
    Loop(j,0,4) {iend[j]=readint(fpin); iend[j]=max(1,min(i4dp[j],iend[j])); }
    printout("normal","ie[] %d %d %d %d, ",iend[0],iend[1],iend[2],iend[3]);
    Loop(j,0,4) {ifix[j]=readint(fpin); ifix[j]=max(0,min(i4dp[j]-1,ifix[j]-1)); }
    printout("normal","ifix[] %d %d %d %d, ",ifix[0],ifix[1],ifix[2],ifix[3]);
    flowoutprop[i]=readdouble(fpin);
    printout("normal","value %lg\n",flowoutprop[i]);
    /* check for consistency */
    L=*dL; if (L<0) L=-L;     /* check dL */
    if (iend[L-1]-ist[L-1]!=1) 
    {
      printout("error exitinit","error, need to specify ist(%d)=iend(%d)\n",L,L); 
      exitm4d(0); 
    }
    if (ist[L-1]==0 && *dL!=L) {*dL=L; printout("normal","dL corrected to %d\n",*dL); }
    if (ist[L-1]==i4d[L-1] && *dL!=-L) {*dL=-L; printout("normal","dL corrected to %d\n",*dL);}
    ipt=In4(ifix,i4dp);  /*check ifix */
    if (matchpc[ipt]>0)
    {
      ipt=matchpc[ipt];
      if (matchpc[ipt]==-1) 
      { 
        iexpand(ipt,i4dp,ii);
        Loop(j,0,3) ifix[j]=ii[j];
        printout("normal c_exitinit","ifix corrected to %d %d %d\n",ifix[0]+1,ifix[1]+1,ifix[2]+1);
      }
      else 
      {
        printout("error "," error ifix point is not an independent point\n");
        exitm4d(0);
      }
    }
    
    n+=16;
  }
}
/* ------------------------------ */
/* fixup  left hand sides for exit equations */
/* call after contcpcfixed if used, or contcpcdu if not */
void c_contcpcexit(FILE *fpin, FILE *fprint)
{ 
  int *i4d, *wherep,*noindppts;  /* need */
  int *flowout; double *flowoutprop;  /* need */  int *wherefw, *cpflop_n, **cpflop_i; double *cam,*cplus,*rho, **cpflop_c;  /* flow rate fixed */
  int *itrange; double *roundoff;/* use if available */
  int *cpc_n, **cpc_i;  double **cpc_c;  /* modify */
  
  int *ipeqexit; /* create */
  int *type,*dL,*ist,*iend,*ifix;
  int i,L,sign,ii[4],ipt,iptfix,iptmfix,iregion,nregion,iim[4],iptm,ieq,i4dp[4];
  
  int ig[4],L2,L3,ib,ic,nfpt,ipt4[4],ieqfix,nadd;   /* for fixed flow rate */
  double area[4][3],radcam[3],*cadd,tol=1.e-8,cadd4[4]={1,-1,-1,1},cfadd4[4]; 
  int nopeq=0,*ipeqtemp;
  
  flowout=(int *)find("flowout");
  if (flowout==0) return;
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  noindppts=(int*)need("noindppts");
  itrange=(int *)find("itrange");
  wherep=(int *)need("wherep"); 
  flowoutprop=(double *)need("flowoutprop");
  cpc_n=(int*)need("cpc_n"); 
  cpc_i=(int**)need("cpc_i"); 
  cpc_c=(double**)need("cpc_c");
  /* for fixed flow rate */
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  wherefw=(int *)need("wherefw");
  cam=(double *)need("cam");
  cplus=(double *)need("cplus");
  rho=(double *)need("rho");
  cpflop_n=(int*)need("cpflop_n"); 
  cpflop_i=(int**)need("cpflop_i"); 
  cpflop_c=(double**)need("cpflop_c");
  geomcinit(); 
  ipeqtemp=(int*)tmalloca(noindppts[0],'i');
  
  printout("normal"," %d outflow regions:\n",flowout[0]);
  nregion=1;
  Loop(iregion,0,flowout[0])
  {
    type=flowout+nregion; dL=type+1; ist=dL+1; iend=ist+4; ifix=iend+4; 
    printout("normal"," exit %d flowout",iregion);
    Loop(i,0,14) printout("normal"," %d",flowout[nregion+i]); printout("normal","\n");
    L=*dL-1; sign=1; if (*dL<0) L=(-*dL)-1; sign=-1;
    Loop(ii[3],ist[3],iend[3])
    { 
      if (itrange!=0 && (ii[3]<itrange[0] || ii[3]>itrange[1]) )continue;
      Loop(i,0,3) ii[i]=ifix[i];  /* put ifix at current time */
      iptfix=In4(ii,i4dp); 
      ipeqtemp[nopeq]=wherep[iptfix];
      nopeq++;
      ii[L]+=sign;
      iptmfix=In4(ii,i4dp);
      Loop(ii[2],ist[2],iend[2]) Loop(ii[1],ist[1],iend[1]) Loop(ii[0],ist[0],iend[0]) 
      {
        ipt=In4(ii,i4dp);
        ieq=wherep[ipt];
        if ( ieq<0) continue;
        if (cpc_n[ieq]>0) {free(cpc_i[ieq]); free(cpc_c[ieq]); cpc_n[ieq]=0; } /* clean */
        if (ipt==iptfix) continue; /* done below */
        ipeqtemp[nopeq]=ieq;
        nopeq++;
        /* other equations */
        if (type[0]==4) /* dp uniforn */
        { 
          cpc_n[ieq]=2;
          cpc_i[ieq]=(int *)smalloca(2,'i');
          cpc_c[ieq]=(double *)smalloca(2*2,'d');
          Loop(i,0,4) iim[i]=ii[i]; iim[L]+=sign; iptm=In4(iim,i4dp);
          cpc_i[ieq][0]=ipt; 
          cpc_i[ieq][1]=iptfix; 
          cpc_c[ieq][0]=1; cpc_c[ieq][1]=-1; 
          Loop(i,0,2) cpc_c[ieq][i+2]=cpc_c[ieq][i];
        }
        else /*  dpdxi uniforn */
        {
          cpc_n[ieq]=4;
          cpc_i[ieq]=(int *)smalloca(4,'i');
          cpc_c[ieq]=(double *)smalloca(4*2,'d');
          Loop(i,0,4) iim[i]=ii[i]; iim[L]+=sign; iptm=In4(iim,i4dp);
          cpc_i[ieq][0]=ipt; cpc_i[ieq][2]=iptm; 
          cpc_i[ieq][1]=iptfix; cpc_i[ieq][3]=iptmfix;
          cpc_c[ieq][0]=1; cpc_c[ieq][1]=-1; cpc_c[ieq][2]=-1; cpc_c[ieq][3]=1;
          Loop(i,0,4) cpc_c[ieq][i+4]=cpc_c[ieq][i];
        }
      }
      /*  mass flow rate fixed equation */
      
      if (type[0]==3 || type[0]==4)
      { 
        ig[3]=ii[3]; L2=(L+1)%3, L3=(L+2)%3;
        ig[L]=ist[L]-1; if (sign==1) ig[L]=ii[L];
        ieqfix=wherep[iptfix];
        Loop(ii[2],ist[2],iend[2]) Loop(ii[1],ist[1],iend[1]) Loop(ii[0],ist[0],iend[0]) 
        { 
          if (wherep[In4(ii,i4dp)]<0) continue; 
          ig[L2]=ii[L2]-1; if (ig[L2]<0 || ig[L2]>i4d[L2]-2) continue;
          ig[L3]=ii[L3]-1; if (ig[L3]<0 || ig[L3]>i4d[L3]-2) continue;
          geomcavector(area,ig,L);
          Loop(ib,0,2) Loop(ic,0,2)
          { 
            ig[L2]=ii[L2]-1+ib; ig[L3]=ii[L3]-1+ic;
            ipt=In4(ig,i4d);
            nfpt=wherefw[ipt];
            if (nfpt<0) continue;
            if (cplus[nfpt]<=0) continue;
            Loop(i,0,3) radcam[i]=rho[ipt]*area[ib+2*ic][i]*sign/(cam[nfpt]+cplus[nfpt]);
            nadd=cpflop_n[nfpt];
            if (nadd==0) continue;
            cadd=(double*)tmalloca(nadd,'d');
            Loop(i,0,nadd) cadd[i]=-(radcam[0]*cpflop_c[nfpt][i]
                                     +radcam[1]*cpflop_c[nfpt][i+nadd]+radcam[2]*cpflop_c[nfpt][i+2*nadd]);
            coefcadd(ieqfix,cpc_n,cpc_i,cpc_c,2,0,nadd,cpflop_i[nfpt],cadd);
          }
        }
        /* coefprint1(ieqfix,cpc_n,cpc_i,cpc_c,2,fprint); */ /* check for now */
        coefcombine(ieqfix,cpc_n,cpc_i,cpc_c,2,tol);
        /* coefprint1(ieqfix,cpc_n,cpc_i,cpc_c,2,fprint); */ /* check for now */
        ipt4[1]=iptfix; ipt4[3]=iptmfix;
        ib=cpc_n[ieqfix];
        Loop(i,0,ib) /* eliminate all but iptfix at exit plane */
        { 
          if (cpc_i[ieqfix][i]==iptfix) continue;
          iexpand(cpc_i[ieqfix][i],i4dp,ii);
          if (ii[L]!=ist[L]) continue;
          ipt4[0]=cpc_i[ieqfix][i];
          ii[L]+=sign;
          ipt4[2]=In4(ii,i4dp);
          Loop(ic,0,4) cfadd4[ic]= -cadd4[ic]*cpc_c[ieqfix][i];
          /*printout("normal","ipt4 %d %d %d %d  cfadd4 %lg %lg %lg %lg\n",ipt4[0],ipt4[1],ipt4[2],ipt4[3],
           cfadd4[0],cfadd4[1],cfadd4[2],cfadd4[3]); */
          if (type[0]==3) coefcadd(ieqfix,cpc_n,cpc_i,cpc_c,2,0,4,ipt4,cfadd4);
          else if (type[0]==4) coefcadd(ieqfix,cpc_n,cpc_i,cpc_c,2,0,2,ipt4,cfadd4);
        }
        coefcombine(ieqfix,cpc_n,cpc_i,cpc_c,2,tol);
        coefcenterfirst(iptfix,ieqfix,cpc_n,cpc_i,cpc_c,2);
        /* coefprint1(ieqfix,cpc_n,cpc_i,cpc_c,2,fprint); */ /* check for now */
      }
    } /* ii[3] */
    nregion+=16;
  }
  if (nopeq>0) 
  { ipeqexit=(int*)createarray("ipeqexit",nopeq+1,'i',0);
    ipeqexit[0]=nopeq;
    Loop(i,0,nopeq) ipeqexit[i+1]=ipeqtemp[i];
  }
}
/* ------------------------------ */
/* fixup right hand sides for exit equations */
/* call after contrhsu */
void c_contrhsexit(FILE *fpin, FILE *fprint)
{
  int *i4d, *wherep;  /* need */
  int *flowout; double *flowoutprop, *pp;  /* need */
  double *rho,*u[3]; /* need for flow rate exit */
  int *itrange; /* use if available */
  double *rhsc;  /* modify */
  
  int *type,*dL,*ist,*iend,*ifix;
  
  int i,L,sign,ii[4],ipt,iptfix,iptmfix,iregion,nregion,iim[4],iptm,ieq,i4dp[4];
  int ieqfix,ig[4],L2,L3,ib,ic;  /* for flow rate exit */
  double area[4][3];
  
  flowout=(int *)find("flowout");
  if (flowout==0) return;
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  itrange=(int *)find("itrange");
  wherep=(int *)need("wherep"); 
  flowoutprop=(double *)need("flowoutprop");
  rhsc=(double *)need("rhsc"); 
  pp=(double *)need("pp");
  
  rho=(double *)need("rho"); /* for flow rate exit */
  u[0]=(double *)need("U1"); 
  u[1]=(double *)need("U2"); 
  u[2]=(double *)need("U3"); 
  
  printout("normal"," %d outflow regions:\n",flowout[0]);
  nregion=1;
  Loop(iregion,0,flowout[0])
  { 
    type=flowout+nregion; dL=type+1; ist=dL+1; iend=ist+4; ifix=iend+4; 
    printout("normal"," exit %d flowout",iregion);
    Loop(i,0,14) printout("normal"," %d",flowout[nregion+i]); printout("normal","\n");
    L=*dL-1; sign=1; if (*dL<0) L=(-*dL)-1; sign=-1;
    Loop(ii[3],ist[3],iend[3])
    {
      if (itrange!=0 && (ii[3]<itrange[0] || ii[3]>itrange[1]) )continue;
      Loop(i,0,3) ii[i]=ifix[i];  /* put ifix at current time */
      iptfix=In4(ii,i4dp);
      ieqfix=wherep[iptfix];
      rhsc[ieqfix]=0;
      if (type[0]==3 || type[0]==4) rhsc[ieqfix]=-flowoutprop[iregion];
      ii[L]+=sign;
      iptmfix=In4(ii,i4dp);
      ig[3]=ii[3]; 
      L2=(L+1)%3, L3=(L+2)%3;
      ig[L]=ist[L]-1; if (sign==1) ig[L]=ii[L];
      Loop(ii[2],ist[2],iend[2]) Loop(ii[1],ist[1],iend[1]) Loop(ii[0],ist[0],iend[0]) 
      { 
        ipt=In4(ii,i4dp);
        ieq=wherep[ipt];
        if ( ieq<0) continue;
        if (ipt!=iptfix) /* type 1 for now */
        {
          Loop(i,0,4) iim[i]=ii[i]; iim[L]+=sign; iptm=In4(iim,i4dp);
          rhsc[ieq]=-pp[ipt]+pp[iptm]+pp[iptfix]-pp[iptmfix];
        }
        if (type[0]==3 || type[0]==4) /* flow rate */
        {
          ig[L2]=ii[L2]-1; ig[L3]=ii[L3]-1;
          geomcavector(area,ig,L);
          Loop(ib,0,2) Loop(ic,0,2)
          { 
            ig[L2]=ii[L2]-1+ib; ig[L3]=ii[L3]-1+ic;
            ipt=In4(ig,i4d);
            rhsc[ieqfix]-=sign*(area[ib+2*ic][0]*rho[ipt]*u[0][ipt]
                                +area[ib+2*ic][1]*rho[ipt]*u[1][ipt]
                                +area[ib+2*ic][2]*rho[ipt]*u[2][ipt]);
          }
        }
      }
      if (type[0]==3 || type[0]==4)
        printout("normal","rhsc at ifix %lg flow rate %lg\n",rhsc[ieqfix],flowoutprop[iregion]);
    } /* ii[3] */
    nregion+=16;
  }
}



