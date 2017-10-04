/* contains tdmaalloc, tdmafree, tdmafromcoef, tdmasolve, tdmaijkcoef, tdmar*/
#include "global.h" 
#include "tdmasubs.h"
#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* real definition in tdmasubs.h */
#if 0
/* tdma structure  commented out here, official definition in tdmasubs.h */ 
/* use cplus instead of cc for coef of ic */
typedef struct 
{ 
  int icmax; /* number of variables */
  int *ic;   /* relative equation indices of variables */
  double *cm; /* coefficient of ic-1 */
  double *cp; /* coefficient of ic+1 */
  void *next; /* address of next Tdma */
}  Tdma;
#endif

/* ----------------------- */
Tdma* tdmaalloc(int i)
{
  Tdma *td; 
  td=(Tdma *)smalloc(sizeof(Tdma));
  td->icmax=i;
  td->ic=(int*)smalloc(i*sizeof(int));
  td->cm=(double*)smalloc(i*sizeof(double));
  td->cp=(double*)smalloc(i*sizeof(double));
  td->next=0;
  return (td);
}
/* ------------------------*/
void tdmafree(Tdma *td)
{
  Tdma *tdn;
  while (td!=0)
  { 
    tdn=(Tdma *)td->next;
    free(td->ic); 
    free(td->cm); 
    free(td->cp);
    free(td);
    td=tdn;
  }
}
/* -------------------------- */
Tdma *tdmafromcoef(int ns, int ne, int *coefn, int **coefi, double **coefc, 
                   double *cplus, int *where)
{ 
  int i,j,imin,icmin=0,*icp, *icm, icmax;    /* to organize tdma strings */
  int nyset,n;
  double *cc, *cp, *cm, cmin;
  int debugb=0;
  
  Tdma *td=0, *td1,*tdnew;                    /* tdma */
  int *iset,is,ie,jtot,jj,itdtot=0,*itry;
  
  n=ne-ns;
  printout("normal","      tdmafromcoef: %d equations,",n);
  /*------------------------organize into strings----------*/  
  icp=(int*)smalloc(n*sizeof(int));
  icm=(int*)smalloc(n*sizeof(int));
  cm=(double*)smalloc(n*sizeof(double));
  cp=(double*)smalloc(n*sizeof(double));
  cc=(double*)smalloc(n*sizeof(double));
  iset=(int*)smalloc(n*sizeof(int));
  /*-----------------------set icp icm cm cp cc---------------*/
  for (i=0;i<n;i++) { icp[i]=-1; icm[i]=-1; cc[i]=0; cm[i]=0; cp[i]=0;}
  Loop(i,0,n)
  {
    icmax=coefn[i+ns];
    if (icmax<=0 || cplus[i+ns]<=0) cc[i]=0;
    else if (icmax==1) cc[i]=cplus[i+ns];     
    else
    { 
      int *itry;
      itry=(int*)smalloc(icmax*sizeof(int));
      cmin=1; icmin=1;
      Loop(j,1,icmax)
      {
        itry[j]=0;
        if (coefc[i+ns][j]>0) continue;
        else 
        {
          if (coefc[i+ns][j]<cmin) {cmin=coefc[i+ns][j]; icmin=j;}
        }
      } /* j */ 
      cc[i]=cplus[i+ns];
      nyset=0;
      while (cmin<0)
      {
        imin=where[coefi[i+ns][icmin]]-ns;
        itry[icmin]=1;
        if (imin>=0 && imin<n && icp[imin]!=i && coefn[imin+ns]>0)
        {
          if (icm[imin]<0) {nyset=1; break;} /* neighbor not set, ok to set */
          if (cmin/cc[i] < cp[icm[imin]]/cc[icm[imin]])
          {nyset=1; icp[icm[imin]]=-1; break;} /* break other chain */
        }
        cmin=1;  /* find next lowest coefficient */
        Loop(j,0,icmax)
        if (coefc[i+ns][j]<cmin && itry[j]==0) {cmin=coefc[i+ns][j]; icmin=j;}
      }
      if (nyset==1)     /* set p side to point to imin */
      {
        cp[i]=cmin; icp[i]=imin;
        icm[imin]=i; cm[imin]=0.;
        Loop(j,1,coefn[imin+ns])
        if (coefi[imin+ns][j]==coefi[i+ns][0]) 
        {cm[imin]=coefc[imin+ns][j]; break;}
        if(debugb>4)     printout("normal","set %d to %d\n",i,imin);
      }
      free(itry);
    } /* icmax>1 */
  } /* i */
  /*------------------------try again-----------------------*/
  Loop(i,0,n)
  { 
    icmax=coefn[i+ns];
    if (icmax>1 && icp[i]==-1)
    { 
      itry=(int*)malloc(icmax*sizeof(int));
      cmin=1;
      Loop(j,0,icmax)
      { 
        itry[j]=0;
        if (coefc[i+ns][j]<cmin) {cmin=coefc[i+ns][j]; icmin=j;}
      } /* j */ 
      nyset=0;
      while (cmin<0)
      {
        imin=where[coefi[i+ns][icmin]]-ns;
        itry[icmin]=1;
        if (imin>=0 && imin<n && icp[imin]!=i && coefn[imin+ns]>0)
        {
          if (icm[imin]<0) {nyset=1; break;} /* neighbor not set, ok to set */
          if (cmin/cc[i] < cp[icm[imin]]/cc[icm[imin]])
          {nyset=1; icp[icm[imin]]=0; break;} /* break other chain */
        }
        cmin=1;  /* find next lowest coefficient */
        Loop(j,0,icmax)
        if (coefc[i+ns][j]<cmin && itry[j]==0) {cmin=coefc[i+ns][j]; icmin=j;}
      }
      if (nyset==1)     /* set p side to point to imin */
      {
        cp[i]=cmin; icp[i]=imin;
        icm[imin]=i; cm[imin]=0.;
        Loop(j,1,coefn[imin+ns])
        if (coefi[imin+ns][j]==coefi[i+ns][0]) {cm[imin]=coefc[imin+ns][j]; break;}
        if(debugb>4)     printout("normal","set %d to %d\n",i,imin);
      }
      free(itry);
    } /* icmax>1 */
  }/* i */
  
  printout("normal"," strings organized,");
  
  /*-----------------------set to tdma structures------------*/
  td1=0;
  for (i=0;i<n;i++) {iset[i]=0;   if (cc[i]==0) iset[i]=1;}  /* omit non-equations */
  for (i=0;i<n;i++) if (iset[i]==0) /* test if already in tdma chain or not wanted*/
  { 
    is=i; jtot=1;     /* find 'is' that starts chain, jtot=chain count */
    while(1) 
    { 
      if (icm[is]<0)             /* found start go to end*/
		{
        ie=i; 
        while (icp[ie]>0) 
        {
          ie=icp[ie]; jtot++; if (debugb>3)  printout("normal"," addp ie jtot %d %d\n",ie,jtot);  
        } 
        break;  /* end while statement */
		}
      if (icm[is]==i)  /* circular loop found */
      {
        ie=i; if  (debugb>3)     printout("normal","circular loop \n"); 
        break; /* end while statement */
      }        
      is=icm[is];
      jtot++; if (debugb>3)  printout("normal"," addm is jtot %d %d\n",is,jtot);
    }
    tdnew=tdmaalloc(jtot);      /* set new tdma string */ 
    tdnew->ic[0]=is;
    /*tdnew->cc[0]=cc[is]; */
    tdnew->cp[0]=cp[is];
    tdnew->cm[0]=0;
    iset[is]=1;
    for (j=1;j<jtot;j++)
    { 
      jj=icp[tdnew->ic[j-1]];
      tdnew->ic[j]=jj;
      tdnew->cm[j]=cm[jj];
      /*tdnew->cc[j]=cc[jj]; */
      tdnew->cp[j]=cp[jj];
      iset[jj]=1;
    }
    tdnew->cp[jtot-1]=0;
    if (debugb>3) printout("normal","td1 td tdnew %lu %lu %lu\n",
                         (unsigned long int)td1,(unsigned long int)td,(unsigned long int)tdnew);
    if (td1==0) {td1=tdnew; td=tdnew; itdtot=0;}
    else {td->next=tdnew; td=tdnew; itdtot++;}
  } /* iset[i]==0 */
  printout("normal"," %d tdma's set\n",itdtot+1);
  
  /* check and print */
  if (debugb>2)
  {
    tdnew=td1; itdtot=0; ie=0;
    while (tdnew!=0)
    {
      itdtot++;
      printout("normal","idtot %d icmax %d\n",itdtot,tdnew->icmax);
      Loop(i,0,tdnew->icmax)
      {
        ie++;
        printout("normal","ic %d cm %lg cc %lg cp %lg coefn %d  ie %d",
               tdnew->ic[i],tdnew->cm[i],cplus[i+ns],tdnew->cp[i],
               coefn[ns+tdnew->ic[i]],ie);
        if (coefn[ns+tdnew->ic[i]]>0)
          printout("normal"," ipt %d\n",coefi[ns+tdnew->ic[i]][0]);
        else printout("warning tdmafromcoef","  error*****\n");
      }
      if (ie>ne) {printout("error tdmafromcoef","too many points\n"); exitm4d(0);}
      tdnew=tdnew->next;
    }
  }
  free(icp); free(icm); free(cm); free(cp); free(cc); free(iset);
  return (td1);
}
/*--------------------------------------------------*/
void tdmasolve(Tdma *td1, int ns, int **coefi, double *cplus, 
               double *bb, double *x)
{ 
  int max, j, *ic;
  double *cpn, *bbn, *dx, *ccn, *cmn, d;
  Tdma *td;   int jtot=0;
  
  td=td1;
  while (td!=0)
  {
    jtot++;
    max=td->icmax; ic=td->ic;
    if (max==1) 
    {
      if (cplus[ic[0]+ns]>0) 
        x[coefi[ns+ic[0]][0]]+=bb[ic[0]+ns]/cplus[ic[0]+ns];
    }
    else 
    {
      dx=(double*)smalloc(max*5*sizeof(double));
      cpn=dx+max; bbn=cpn+max; ccn=bbn+max; cmn=ccn+max;
      Loop(j,0,max)
      { 
        if (cplus[td->ic[j]+ns]>0) 
        {
          ccn[j]=cplus[ic[j]+ns];
          cmn[j]=td->cm[j];
          cpn[j]=td->cp[j];
          bbn[j]=bb[ic[j]+ns];
        }
        else { ccn[j]=1; cmn[j]=0; cpn[j]=0; bbn[j]=0;}
      }
      if (cmn[0]!=0 || cpn[max-1]!=0)  /* use repeat tdma */
        tdmar(cmn,ccn,cpn,bbn,dx,max);
      else /* reg tdma no repeat */ 
      {
        if (max==2) 
        {
          dx[1] = (bbn[1]-cmn[1]*bbn[0]/ccn[0])/(ccn[1]-cmn[1]*cpn[0]/ccn[0]);
          dx[0] = bbn[0]/ccn[0] -dx[1]*cpn[0]/ccn[0];
        }
        else   /* max>2 */
        {
          cpn[0] = -cpn[0]/ccn[0];
          bbn[0] = bbn[0]/ccn[0];
          for (j=1;j<max-1;j++)
          { 
            d = ccn[j]+cmn[j]*cpn[j-1];
            bbn[j] = (bbn[j]-cmn[j]*bbn[j-1])/d;
            cpn[j] = -cpn[j]/d;
          }
          dx[max-1] = (bbn[max-1]-cmn[max-1]*bbn[max-2])
          /(ccn[max-1]+cmn[max-1]*cpn[max-2]);
          for (j=max-2;j>=0;j--) dx[j] = bbn[j]+cpn[j]*dx[j+1];
        }
      }
     
      Loop(j,0,max) x[coefi[ns+ic[j]][0]]+=dx[j];
      free(dx);
    }
    td=td->next;
  }
}

/*--------------------------------------------------*/
Tdma *tdmaijkcoef(int L, int ns, int ne, int *i4d, int *coefn, int **coefi, 
                  double **coefc, double *cplus, int *where, int *match)
{	
  int i,j,ii[4],idmax,n,icmax,ipt,i1,ieq,itdtot=0,noeqstot=0;
  int *ic,*idone;
  Tdma *td1=0,*td=0,*tdnew;
  
  idmax=i4d[L];
  ic=(int*)smalloc(idmax*sizeof(int));
  idone=(int*)smalloc((ne-ns)*sizeof(int));
  Loop(i,0,ne-ns) idone[i]=0;
  
  Loop(n,ns,ne) /* look for starts */
  if (idone[n-ns]==0)
  { 
    icmax=0;
    if (coefn[n]<=0 || cplus[n]<=0) continue;
    ipt=coefi[n][0];
    /*printout("normal","start eq %d pt %d ",n,ipt);*/
    ic[0]=n;
    idone[n-ns]++;
    iexpand(ipt,i4d,ii);
    i1=ii[L]+1;
    /*printout("normal","ii %d %d %d %d L %d i1 %d idmax %d more:\n",ii[0],ii[1],ii[2],ii[3],L,i1,idmax);*/
    Loop (ii[L],i1,idmax) /* find eq no for next valid point */
    {  
      ipt=In4(ii,i4d);
      if (match[ipt]>=0) ipt=match[ipt];
      ieq=where[ipt];
      if (ieq<0) continue;
      if (ieq==ic[icmax]) continue;	
      if (cplus[ieq]<=0) continue;
      if (coefn[ieq]<=0) continue;
      if (idone[ieq-ns]>0) continue;
      icmax++;
      ic[icmax]=ieq;
      idone[ieq-ns]++;
      /*printout("normal"," %d %d %d,",ipt,ieq,icmax);*/
    }
    /*printout("normal","\n");*/
    icmax++;   /* ok have equation list */
    tdnew=tdmaalloc(icmax);
    tdnew->icmax=icmax;
    Loop(i,0,icmax)
    {
      ieq=ic[i];
      noeqstot++;
      tdnew->ic[i]=ieq;
      tdnew->cm[i]=0;
      tdnew->cp[i]=0;
      if (i>0)
        Loop(j,0,coefn[ieq])
        if (where[coefi[ieq][j]]==ic[i-1]) 
        {tdnew->cm[i]=coefc[ieq][j]; break;}
      if (i+1<icmax)
        Loop(j,0,coefn[ieq])
        if (where[coefi[ieq][j]]==ic[i+1]) 
        {tdnew->cp[i]=coefc[ieq][j]; break;}
    }
    /* check for repeat */
    ieq=ic[0];
    Loop(j,0,coefn[ieq])
    if (where[coefi[ieq][j]]==ic[icmax-1])
    {tdnew->cm[0]=coefc[ieq][j]; break;}
    ieq=ic[icmax-1];
    Loop(j,0,coefn[ieq])
    if (where[coefi[ieq][j]]==ic[0])
    {tdnew->cp[icmax-1]=coefc[ieq][j]; break;}
    
    if (td1==0) {td1=tdnew; td=tdnew;}
    else {td->next=tdnew; td=tdnew;}
    itdtot++;
    /*printout("normal","itdtot %d noeqstot %d\n",itdtot,noeqstot);*/
  }  
  free(idone); free(ic); 
  printout("normal"," L=%d, %d strings for %d equations\n",L,itdtot,noeqstot);
  return td1;
}
/* --------------------------------- */
/* tdma for repeating points, 
 eqs solved:  cm[j]*x[j-1]+cc[j]*x[j]+cp[j]*x[j+1]=bb[j], j=0,,j<max
 cm of first point is coef of last pt in first eq. 
 cp of last pt is coef of first pt in last eq.
 note the incoming arrays are destroyed, if needed copy before calling 
 solution is returned in x
 */
void tdmar(double *cm, double *cc, double *cp, double *bb, double *x, int max)
{  
  double cmr,ccr,cpr,bbr,d;
  int j;
  
  cmr=cm[0]; ccr=cc[0]; cpr=cp[0]; bbr=bb[0];
  bb[0]=0; cp[0]=0; cc[0]=1;
  
  for (j=1;j<max;j++)
  {
    d=cc[j]+cm[j]*cp[j-1];
    /* if (d<.001*cm[j]) d=.001*cm[j];*/   /* safety check */
    bb[j]=(bb[j]-cm[j]*bb[j-1])/d;
    cp[j]=-cp[j]/d;
    cc[j]=-(cm[j]*cc[j-1])/d;
  }
  
  cc[max-1]=cc[max-1]+cp[max-1];
  x[max-1]=bb[max-1];
  for (j=max-2;j>0;j--)
  {
    cc[j]=cc[j]+cp[j]*cc[j+1];
    x[j]=bb[j]+cp[j]*x[j+1];
  }
  
  x[0]=(bbr-cmr*x[max-1]-cpr*x[1])/(ccr+cmr*cc[max-1]+cpr*cc[1]);
  for (j=1;j<max;j++)
  {
    x[j]=x[j]+cc[j]*x[0];
  }
}

