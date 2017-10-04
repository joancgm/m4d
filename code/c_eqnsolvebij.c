/* containd c_eqnsolvebij */
#include "global.h"
#include "tdmasubs.h"
#include "eqnsubs.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

double fixdbij(double *bij, double *dbij, double gtol);
/* solve the set of bij equations of the form c_n, c_i, c_c  neqs,
 where (interconnect info)    with  bij eqn source term =
 +beqm[ij] + beqd[ij][nm] db[nmcen]+camin*dbcen
 read information on type of solve to be used 
 any point which lack coefficients will be left unchanged 
 */

void c_eqnsolvebij(FILE *fpin, FILE *fprint)
{ 
  double *beqm,*beqd, *bij;   /* needed arrays */
  int *i4d,*npts,*neqs,*where,*cn,**ci, *match, *whois; 
  char *clt,*csym;
  double **cc; 
  int *itrange; double *camin; /* use if  available */
  double *dbij; /* set */
  char *ncamin,oldnewch; int nomit;  /*  input parameters */
  char omit[20];
  char *nameproc; int iterproc;   /* read procedure choices until nameproc = "" */
  
  double *cplus, *rhs6, *csave=0, *ddbij; char *nybijeq; /* temporary arrays */
  Tdma *tdma1[4]={0,0,0,0};
  int istdma[4]={0,0,0,0},jtdma; char ctdma[4]={'t','i','j','k'};
  
  int i,j,k,m,ieqstart,ieqend, iter,ii[4],ieq, noproc;
  double err[18]; int imin[6],imax[6],iimin[4],iimax[4];
  double errrmsold,errrms,errmax,cb;
  double b[6],db[6],freal,fmin=1,gtol=.001,ddbm,ddsq,ddbrms=0,ddbmax,ddbrmsold; 
  int noreal=0,nchange;
  
  ncamin=readname(fpin);
  camin=(double *)find(ncamin);
  if (camin>0) printout("normal","using %s, ",ncamin);
  else 
  {  
    camin=(double *)need("cambij");
    printout("normal","using cambij, ");
  }
  oldnewch=read1charname(fpin); /* read input (procedure read later) */
  nomit=readint(fpin);
  printout("normal","oldnewch=%c, omit %d pt types:",oldnewch,nomit);
  Loop(i,0,nomit) 
  {
    omit[i]=read1charname(fpin); 
    printout("normal"," %c",omit[i]); 
  }
  i4d=(int *)need("idim4d");/* get needed arrays */
  beqm=(double *)need("rhsbij");
  beqd=(double *)need("beqd");
  bij=(double *)need("bij");
  npts=(int *)need("nogpts");  
  neqs=(int*)need("noindfwpts");
  where=(int *)need("wherefw");
  whois=(int *)need("whoisfw");
  match=(int *)need("match");
  cn=(int *)need("coef_n");
  ci=(int **)need("coef_i");
  cc=(double **)need("coef_c");
  clt=(char *)need("clt");
  csym=(char *)need("csym");
  itrange=(int *)find("itrange");
  ieqstart=0; ieqend=neqs[0];
  if (itrange>0) {ieqstart=neqs[itrange[0]+1]; ieqend=neqs[itrange[1]+2]; }
  istdma[0]=ieqstart;
  printout("normal",",  eq range %d, %d\n",ieqstart,ieqend);
  
  dbij=(double *)find("dbij");
  if (dbij==0) 
  { 
    dbij=(double *)createarray("dbij",npts[0]*6,'d',0);
    Loop(i,0,npts[0]*6) dbij[i]=0;
  }
  else if (oldnewch=='n') Loop(i,0,npts[0]*6) dbij[i]=0;
  
  /* n y solve equations  this could be a separate routine to be done once */
  nybijeq=(char *)tmalloca(neqs[0]*4,'c'); 
  Loop(i,0,neqs[0]*4) nybijeq[i]='y';
  Loop(ii[3],0,i4d[3]) Loop(ii[2],0,i4d[2]) Loop(ii[1],0,i4d[1]) Loop(ii[0],0,i4d[0])
  {
    m=In4(ii,i4d);
    ieq=where[match[m]];
    if (ieq<0) continue;
    Loop(j,0,nomit)
    if (clt[m]==omit[j])
      Loop(i,0,4) nybijeq[ieq+i*neqs[0]]='n';
    Loop(j,0,6) if (csym[j]!='n')
    { 
      if (j==0 && ii[0]!=0) continue;
      if (j==1 && ii[0]!=i4d[0]-1) continue;
      if (j==2 && ii[1]!=0) continue;
      if (j==3 && ii[1]!=i4d[1]-1) continue;
      if (j==4 && ii[2]!=0) continue;
      if (j==5 && ii[2]!=i4d[2]-1) continue; 
      if (csym[j]=='x') {nybijeq[ieq+neqs[0]]='n'; nybijeq[ieq+2*neqs[0]]='n';} 
      else if (csym[j]=='y') {nybijeq[ieq+neqs[0]]='n'; nybijeq[ieq+3*neqs[0]]='n';} 
      else if (csym[j]=='z') {nybijeq[ieq+2*neqs[0]]='n'; nybijeq[ieq+3*neqs[0]]='n';}
      /*printout("normal","csym %c j %d ii %d %d %d %d ieq %d nybijeq %c %c %c %c\n",
       csym[j],j,ii[0],ii[1],ii[2],ii[3],ieq,nybijeq[ieq],nybijeq[ieq+neqs[0]],
       nybijeq[ieq+2*neqs[0]],nybijeq[ieq+3*neqs[0]]);
       */
    }
  }
  printout("normal","nybijeq set\n"); 
  /* alter beqd putting into eqs 22 33 12 13 23  db11=-db22-db33 */
  Loop(i,ieqstart,ieqend)
  { 
    Loop(j,1,6)
    cb=beqd[i+j*neqs[0]];
    Loop(k,0,3) beqd[i+(j+6*k)]-=cb;
  }
  
  csave=(double *)tmalloca(neqs[0],'d');	
  Loop(i,ieqstart,ieqend) /* save center coefficient and add camin */
  if (cn[i]>0)  
  { 
    csave[i]=cc[i][0]; 
    if (camin>0) cc[i][0]+=camin[i];
  }
  
  /* initial calcs and analysis */
  cplus=(double *)tmalloca(neqs[0],'d');	 
  eqncplus(ieqstart,ieqend,cn,cc, cplus); /* cplus=sum of + coef (space) */
  printout("normal","cplus set \n"); 
  rhs6=(double *)tmalloca(neqs[0]*6,'d'); /* initial residual */
  eqnrhsbij(ieqstart,ieqend,cn,ci,cc,npts[0],neqs[0],
            beqm,dbij,beqd,nybijeq,rhs6); 
  
  printout("normal"," rhs6 set \n"); 
  errrmsold=0; errmax=0;  /* initial error analysis usint cplus */
  Loop(j,0,6) 
  { 
    eqnerrorwho(ieqstart,ieqend,rhs6+j*neqs[0],cplus,nybijeq+max(0,(j-2)*neqs[0]),err+3*j,imin+j,imax+j);
    iexpand(whois[imin[j]],i4d,iimin); iexpand(whois[imax[j]],i4d,iimax);
    printout("normal"," bij %d:  err-rms %lg    err-min %lg at %d %d %d %d   err-max %lg at %d %d %d %d\n",
            j,err[3*j],err[3*j+1],iimin[0],iimin[1],iimin[2],iimin[3],
            err[3*j+2],iimax[0],iimax[1],iimax[2],iimax[3]);
    errrmsold+=err[3*j]*err[3*j];
    errmax=max(errmax,max(-err[3*j+1],err[3*j+2]));
  }
  errrmsold=sqrt(errrmsold/6.);
  printout("normal"," init:  err-rms %lg    |err-max| %lg\n",errrmsold,errmax);
  
  ddbij=(double *)tmalloca(neqs[0]*6,'d'); /* for change analysis */
  Loop(i,0,neqs[0]*6) ddbij[i]=0;
  
  /*----------------------------  read and do choices --------------------*/
  while(1)
  { 
    nameproc=readname(fpin);
    noproc=strlen(nameproc);
    printout("normal","Procedure: %s",nameproc);
    if (nameproc[0]=='e') { printout("normal","\n"); break; }  /* end */
    iterproc=readint(fpin);
    printout("normal"," ,  %d iterations\n",iterproc);
    if (iterproc==0) break;    /* end */
    
    Loop(i,0,noproc) /* check for need of tdma's */
    {	
      if (nameproc[i]==ctdma[0] && tdma1[0]==0) 
        tdma1[0]=tdmafromcoef(ieqstart,ieqend,cn,ci,cc,cplus,where);
      Loop(k,1,4)
      if (nameproc[i]==ctdma[k] && tdma1[k]==0) 
        tdma1[k]=tdmaijkcoef(k-1,ieqstart,ieqend,i4d,cn,ci,cc,cplus,where,match);
    }
    Loop(iter,0,iterproc)   /* do solves */
    { 
      Loop(ieq,ieqstart,ieqend)  /* save dbij for ddb analysis */
      Loop(j,0,6) ddbij[ieq+j*neqs[0]]=-dbij[whois[ieq]+j*npts[0]];
      
      Loop(i,0,noproc)
      {  				
        if (nameproc[i]=='c') 
          eqncenterbij(ieqstart,ieqend,neqs[0],npts[0],cplus,beqd,rhs6,ci,nybijeq,bij,dbij);
        else 
        {  
          jtdma=-1; Loop(k,0,4) if (nameproc[i]==ctdma[k]) jtdma=k;
          if (jtdma==-1) 
          { 
            printout("normal"," procedure %c not found\n",nameproc[i]); 
            if (noproc==1) break; 
          }
          else Loop(j,1,6) 
          { 
            /* fix up cplus for sym planes */
            Loop(ieq,ieqstart,ieqend) 
            if (nybijeq[ieq+max(0,j-2)*neqs[0]]=='n') cplus[ieq]=-cplus[ieq];
            tdmasolve(tdma1[jtdma],istdma[jtdma],ci,cplus,rhs6+j*neqs[0],dbij+j*npts[0]);
            /* fix up cplus for sym planes */
            Loop(ieq,ieqstart,ieqend) 
            if (nybijeq[ieq+max(0,j-2)*neqs[0]]=='n') cplus[ieq]=-cplus[ieq];
          }
          Loop(ieq,ieqstart,ieqend) /* set db11 */
          { 
            j=whois[ieq];
            dbij[j]=-dbij[j+npts[0]]-dbij[j+2*npts[0]]-bij[j]
            -bij[j+npts[0]]-bij[j+2*npts[0]];
          }
        }
        /* fix dbij for realizability */
        fmin=1; noreal=0;
        Loop(ieq,ieqstart,ieqend)
        { 
          Loop(j,0,6) {b[j]=bij[whois[ieq]+npts[0]*j]; db[j]=dbij[whois[ieq]+npts[0]*j]; } 
          freal=fixdbij(b,db,gtol);
          if (freal<1) 
          { 
            Loop(j,0,6) dbij[whois[ieq]+npts[0]*j]=db[j];
            if (freal<fmin) fmin=freal;
            noreal++;
          }
        }
        eqnrhsbij(ieqstart,ieqend,cn,ci,cc,npts[0],neqs[0],
                  beqm,dbij,beqd,nybijeq,rhs6); 
      }
      ddbrmsold=ddbrms;  /* ddbij analysis */
      ddbrms=0; ddbmax=0; nchange=0; 
      Loop(ieq,ieqstart,ieqend)  
      {  
        ddbm=0; ddsq=0;
        Loop(j,0,6) 
        { 
          ddbij[ieq+j*neqs[0]]+=dbij[whois[ieq]+j*npts[0]];
          ddbm=max(ddbm,ddbij[ieq+j*neqs[0]]);
          ddsq+=ddbij[ieq+j*neqs[0]]*ddbij[ieq+j*neqs[0]];
        }
        if (ddbm>0)
        {  
          nchange++;
          ddbmax=max(ddbm,ddbmax);
          ddbrms+=ddsq;
        }
      }
      if (iter==0 || iter==iterproc-1 || ddbrms>ddbrmsold)
        printout("normal","  iter %d  ddb-rms %lg  ddb-max %lg  for %d of %d eqs realfix %d pts fmin %lg\n",
                iter,sqrt(ddbrms/(6*nchange)),ddbmax,nchange,ieqend-ieqstart,noreal,fmin);
      
      errrms=0; errmax=0;  /*  error analysis usint cplus */
      Loop(j,0,6) 
      { 
        eqnerrorwho(ieqstart,ieqend,rhs6+j*neqs[0],cplus,nybijeq+max(0,(j-2)*neqs[0]),
                    err+3*j,imin+j,imax+j); 
        /*printout("normal"," bij %d:  err-rms %lg    err-min %lg at %d   err-max %lg at %d\n",j,
         err[3*j],err[3*j+1],imin[j],err[3*j+2],imax[j]); */
        errrms+=err[3*j]*err[3*j];
        errmax=max(errmax,max(-err[3*j+1],err[3*j+2]));
      }
      errrms=sqrt(errrms/6.);
      if (iter==0 || iter==iterproc-1 || errrms>errrmsold)
        printout("normal","  iter %d  err-rms %lg    err-max %lg\n",iter,errrms,errmax);
      errrmsold=errrms;
    }
  }
  /* end of solve , clean up*/
  printout("normal","end solve\n"); 
  Loop(i,ieqstart,ieqend)
  if (cn[i]>0) cc[i][0]=csave[i];
  Loop (k,0,4) if (tdma1[k]>0) tdmafree(tdma1[k]);
  /* set matched points */
  Loop(i,0,npts[0]) Loop(j,0,6) dbij[i+j*npts[0]]=dbij[match[i]+j*npts[0]];
}


