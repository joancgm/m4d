/* contains c_contcpcfixed */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

/* continuity equations, pressure-correction coefficients  for between the points volumes */
/* fix  up cpc_c 0 to be mixed coefficients for eqs to be solved, 
 and cpc_c 1 to be coefs for permanent P correction to cont equation */

void c_contcpcfixed(FILE *fpin,FILE *fprint)
{ 
  int *i4d,*noindppts,*wherep,*whoisp;    /* need */
  int *cpc_n,**cpc_i; double **cpc_c; /*  modify */
  int *cpc1_n; /* create */
  int *itrange; double *roundoff;
  int ieq,n,i,j,nn,icen,icena,ia,iaa,ieqa,*cinew,na;
  double *ccnew,cadd;
  int neqs,neqe;
  int ipose[11]; 
  double cen,sump,pose=0,disp,posemin,dispmin,facmix,facmax,cmin;
  int db=1; /* change to 1 or 2 for debug prints */
  int iperr[4],i4dp[4],iter;
  double posel[2],displ[2];
  double tol=1.e-8;
  
  double facrc,facrcperm,poselim,displim,poselimperm,displimperm;   /* input parameters */
  
  facrc=readdouble(fpin);
  facrcperm=readdouble(fpin);
  poselim=readdouble(fpin);
  displim=readdouble(fpin);
  poselimperm=readdouble(fpin);
  displimperm=readdouble(fpin);
  facrcperm=min(facrcperm,facrc);
  poselim=min(poselim,1);
  displim=min(displim,-1);
  poselimperm=min(poselim,poselimperm);
  displimperm=min(displim,displimperm);
  
  printout("normal"," fix cpc: facrc %lg, facrcperm %lg, poselim %lg, displim %lg",
          facrc,facrcperm,poselim,displim);
  printout("normal",", poselimperm %lg, displimperm %lg\n",poselimperm,displimperm);
  
  posel[0]=poselimperm; posel[1]=poselim;
  displ[0]=displimperm; displ[1]=displim;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  wherep=(int *)need("wherep");
  whoisp=(int *)need("whoisp");
  noindppts=(int *)need("noindppts");
  cpc_n=(int *)need("cpc_n");
  cpc_i=(int **)need("cpc_i");
  cpc_c=(double **)need("cpc_c");
  nn=noindppts[0];
  itrange=(int *)find("itrange");
  if (itrange==0) {neqs=0; neqe=nn; }
  else { neqs=noindppts[itrange[0]+1];  neqe=noindppts[itrange[1]+2]; }
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  cpc1_n=(int *)createarray("cpc1_n",noindppts[0],'i',0);
  Loop(i,0,noindppts[0]) cpc1_n[i]=cpc_n[i];
  /* check consistency of ci indices */
  j=0;
  Loop(ieq,neqs,neqe)
  { 
    n=cpc_n[ieq];
    if (n>0)    Loop(i,0,n)  if (wherep[cpc_i[ieq][i]]<0) j++;
  }
  if (j>0) 
  {
    printout("error c_contcpcfixed"," %d bad pc indices in cpc_i, stop \n",j); 
    exitm4d(0); 
  }
  /* add in rc   to c 0*/
  Loop(ieq,neqs,neqe)
  { 
    n=cpc_n[ieq];
    if (n>0)
    { 
      Loop(i,0,n) cpc_c[ieq][i]+=facrc*cpc_c[ieq][i+n];
    }
  }
  /* scale c 1 to permanent part */
  Loop(ieq,neqs,neqe)
  { 
    n=cpc_n[ieq];
    if (n>0)
    { 
      Loop(i,0,n) cpc_c[ieq][i+n]*=facrcperm;
    }
  }  
  /* check coefficients sum,pose,disp */
  j=0; Loop(i,0,11) ipose[i]=0; posemin=1; dispmin=0;
  Loop(ieq,neqs,neqe) if (cpc_n[ieq]>1)
  { 	
    cen=cpc_c[ieq][0]; 
    if (cen<=0) {pose=-10; disp=-10; ipose[0]++; }
    else
    {  
      icen=0;
      cen=cpc_c[ieq][icen]; sump=0; cmin=0;
      n=cpc_n[ieq];
      Loop(i,0,n) 
      {  
        if (cpc_c[ieq][i]>0) sump+=cpc_c[ieq][i];
        else cmin=min(cmin,cpc_c[ieq][i]);
      }
      pose=cen/sump;
      disp=cmin/cen;
      i=min(10,pose*10+1); ipose[i]++;
    }
    posemin=min(posemin,pose); dispmin=min(dispmin,disp);
    if (pose<posel[0] || disp<displ[0]) j++;
  }
  printout("normal","perm. fix %d p-eqs, pose min= %lg, disp min=%lg\n",j,posemin,dispmin);
  printout("normal","pose groups");
  Loop(i,0,11) printout("normal","  %d <%lg,",ipose[i],0.1*i);
  printout("normal","\n");
  
  /* look for positive outside coefficients and mix
   iter=0 permanent contribution, iter=1, cpc to be solved */
  Loop(iter,0,2)
  {  
    printout("normal","iter %d poselim %lg displim %lg\n",
            iter,posel[iter],displ[iter]);
    if (posel[iter]<0) continue;
    facmax=0;
    Loop(ieq,neqs,neqe) 
    {  
      n=cpc_n[ieq];
      if (n<2) continue;
      icen=0;
      cen=cpc_c[ieq][icen]; sump=0; cmin=0;
      Loop(i,0,n) 
      {  
        if (cpc_c[ieq][i]>0) sump+=cpc_c[ieq][i];
        else cmin=min(cmin,cpc_c[ieq][i]);
      }
      pose=cen/sump;
      disp=cmin/cen;
      if (pose>.9999) continue;
      if (posel[iter]==1) facmix=1;
      else
      {
        facmix=max(cmin/displ[iter]-cen,sump*posel[iter]-cen)/(sump-cen);
        facmix=max(0,min(1,facmix));
      }
      if (db>1 || facmix>facmax || (db>0 && pose<posel[iter] && disp<displ[iter]) )
      { 
        iexpand(whoisp[ieq],i4dp,iperr);
        printout("normal","ieq %d, %d %d %d %d facmix %lg cen %lg sump %lg cmin %lg pose %lg disp %lg\n",
                ieq,iperr[0],iperr[1],iperr[2],iperr[3],facmix,cen,sump,cmin,pose,disp);
      }
      if (facmix>facmax) facmax=facmix;
      { 
        Loop(ia,0,n) 
        if (cpc_c[ieq][ia]>0 && ia!= icen)
        {
          if (db>2) { Loop(i,0,n) printout("normal"," %lg",cpc_c[ieq][i]); printout("normal","\n"); }
          if (ia<n) 
          {  
            ieqa=wherep[cpc_i[ieq][ia]]; 
            if (db>2) printout("normal","coef %d %d %lg ieqa %d\n",ia,cpc_i[ieq][ia],cpc_c[ieq][ia],ieqa);
            na=cpc_n[ieqa];
            Loop(icena,0,na) if (wherep[cpc_i[ieqa][icena]]==ieqa) break;
            if (icena<na)
            { 
              Loop(iaa,0,na) if (cpc_i[ieqa][iaa]==cpc_i[ieq][icen]) break;
              if (iaa==na) /*  add space in ieqa equation */
              { 
                cinew=(int *)smalloca(na+1,'i');
                ccnew=(double *)smalloca(2*(na+1),'d');
                Loop(i,0,na)
                { 
                  cinew[i]=cpc_i[ieqa][i]; 
                  ccnew[i]=cpc_c[ieqa][i];
                  ccnew[i+na+1]=cpc_c[ieqa][i+na];
                }
                cinew[na]=cpc_i[ieq][icen];
                iaa=na;
                ccnew[na]=0;
                ccnew[na+na+1]=0;
                free(cpc_i[ieqa]); cpc_i[ieqa]=cinew;
                free(cpc_c[ieqa]); cpc_c[ieqa]=ccnew;
                na++; cpc_n[ieqa]=na;
              }
              /* now do additions */
              cadd=facmix*cpc_c[ieq][ia];
              cpc_c[ieq][icen]+=cadd;
              cpc_c[ieqa][icena]+=cadd;
              cpc_c[ieqa][iaa]-=cadd;
              cpc_c[ieq][ia]-=cadd;
              if (iter==0)  /* add to permanent part */
              {	
                cpc_c[ieq][icen+n]+=cadd;
                cpc_c[ieqa][icena+na]+=cadd;
                cpc_c[ieqa][iaa+na]-=cadd;
                cpc_c[ieq][ia+n]-=cadd;
              }
              j++;
            }
          }
        }
      }
    }
    printout("normal","iter %d poselim %lg displim %lg maximum facmix= %lg\n",
            iter,posel[iter],displ[iter],facmax);
  }
  /*  reorder coefs and set cpc1_n, the number of coefs in the first set which are > roundoff */
  Loop(ieq,neqs,neqe)
  if (cpc_n[ieq]>1)
  {   
    cmin=cpc_c[ieq][0]*tol;
    Loop(n,1,cpc_n[ieq]) /* look for a roundoff coef */
    {  
      if (abs(cpc_c[ieq][n]) > cmin) { cpc1_n[ieq]=n+1; continue; }
      Loop(i,n+1,cpc_n[ieq]) /* look for non roundoff one to trade places with */
      {  
        if (abs(cpc_c[ieq][i]) <= cmin) continue;
        j=cpc_i[ieq][i]; cpc_i[ieq][i]=cpc_i[ieq][n]; cpc_i[ieq][n]=j;
        cen=cpc_c[ieq][i]; cpc_c[ieq][i]=cpc_c[ieq][n]; cpc_c[ieq][n]=cen;
        cen=cpc_c[ieq][i+cpc_n[ieq]]; 
        cpc_c[ieq][i+cpc_n[ieq]]=cpc_c[ieq][n+cpc_n[ieq]]; 
        cpc_c[ieq][n+cpc_n[ieq]]=cen;
        cpc1_n[ieq]=n+1;
      }
    }
  }
}
