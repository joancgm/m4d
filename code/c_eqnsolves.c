/* contains c_eqnsolves */
#include "global.h"
#include "tdmasubs.h"
#include "eqnsubs.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

/* solve a set of equations of the form c_n, c_i, c_c  neqs, where (interconnect info)
 read information on type of solve to be used 
 any point which lack coefficients will be left unchanged 
 set up as group of routine which need no input */

void c_eqnsolves(FILE *fpin, FILE *fprint)
{ 
  char *nrhs, *ncam,*nchange;     /*names of  arrays  */
  double  *rhsin, *camin, *ch;  /* corresponding arrays rhs(needed) cam(optional) change */
  
  int nomit;  /* more input parameters */
  char omit[20],oldnewch;
  char *nameproc; int iterproc;   /* read procedure choices until nameproc = "" */
  
  int *i4d,*npts,*neqs,*where,*whois,*cn,**ci, *match; /* other needed arrays  */
  double **cc; char *clt;  
  int *itrange; /* use itrange if  available */
  
  double *rhs,*cplus,*csave=0; /* temporary arrays */
  Tdma *tdma1[4]={0,0,0,0};
  int itdma=4,istdma[4]={0,0,0,0}; 
  char ctdma[4]={'t','i','j','k'};
  
  double err[4]={0,0,0,0},errinit[3];
  int i,j,k,iter,ieqstart,ieqend,itot,itota,noproc;
  char *ny=0;
  int imin,imax,iimin[4],iimax[4];
  
  nrhs=readname(fpin);     /* read input parameters (procedure read later ) */
  ncam=readname(fpin);
  nchange=readname(fpin);
  printout("normal","on the grid equations with rhs=%s, center_add=%s, change=%s\n",nrhs,ncam,nchange);
  oldnewch=read1charname(fpin);
  nomit=readint(fpin);
  printout("normal","oldnew%s %c omit %d pt types:",nchange,oldnewch,nomit);
  Loop(i,0,nomit) {omit[i]=read1charname(fpin); printout("normal"," %c",omit[i]); }
  
  npts=(int *)need("nogpts");   /* get needed arrays */
  i4d=(int *)need("idim4d");
  neqs=(int*)need("noindfwpts");
  where=(int *)need("wherefw");
  whois=(int *)need("whoisfw");
  match=(int *)need("match");
  clt=(char *)need("clt");
  cn=(int *)need("coef_n");
  ci=(int **)need("coef_i");
  cc=(double **)need("coef_c");
  rhsin=(double *)need(nrhs);
  camin=(double *)find(ncam);
  itrange=(int *)find("itrange");
  ieqstart=0; ieqend=neqs[0];
  if (itrange>0) {ieqstart=neqs[itrange[0]+1]; ieqend=neqs[itrange[1]+2]; }
  istdma[3]=ieqstart;
  printout("normal",",  eq range %d, %d\n",ieqstart,ieqend);
  
  ch=(double *)find(nchange);
  if (ch==0 || oldnewch=='n') 
  { 
    ch=(double *)createarray(nchange,npts[0],'d',0);
    Loop(i,0,npts[0]) ch[i]=0;
  }
  /*   end intial array setup */
  /* initial calcs and analysis */
  rhs=(double *)tmalloca(neqs[0],'d');
  cplus=(double *)tmalloca(neqs[0],'d');
  
  if (nomit>0) /* change sign of coef_n so that eqs are omitted */
  { 
    itot=0; itota=0;
    Loop(i,ieqstart,ieqend)
    if (cn[i]>0)
    {
      itota++;
      Loop(j,0,nomit) 
      if (clt[ci[i][0]]==omit[j]) {cn[i]=-cn[i]; itot++; break;}
    }
    printout("normal"," %d of %d active eqs, type deactivated\n",itot,itota);
  }
  if (camin>0)   /* add camin to center coef */
  {
    csave=(double*)tmalloca(neqs[0],'d');
    Loop(i,ieqstart,ieqend)
    if (cn[i]>0)
    {
      csave[i]=cc[i][0];
      cc[i][0]+=camin[i];
    }
  }
  else printout("normal","no center pt additions\n");
  
  eqncplus(ieqstart,ieqend,cn,cc, cplus); /* cplus=sum of + coef */
  /* printout("normal","cplus set \n"); */
  eqnrhs(ieqstart,ieqend,cn,ci,cc,ch,rhsin,  rhs);  /* set rhs with  ch included */
  /* printout("normal"," rhs set \n"); */
  eqnerrorwho(ieqstart,ieqend,rhs,cplus,ny, err,&imin,&imax); /* error analysis */ 
  iexpand(whois[imin],i4d,iimin); iexpand(whois[imax],i4d,iimax);
  printout("normal"," init:  err-rms %lg    err-min %lg at %d %d %d %d    err-max %lg at %d %d %d %d\n",
          err[0],err[1],iimin[0],iimin[1],iimin[2],iimin[3],
          err[2],iimax[0],iimax[1],iimax[2],iimax[3]);
  Loop(i,0,3) errinit[i]=err[i];
  
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
    
    Loop(i,0,noproc)
    {
      if (nameproc[i]==ctdma[0] && tdma1[0]==0) 
        tdma1[0]=tdmafromcoef(ieqstart,ieqend,cn,ci,cc,cplus,where);
      Loop(k,1,itdma)
      if (nameproc[i]==ctdma[k] && tdma1[k]==0) 
        tdma1[k]=tdmaijkcoef(k-1,ieqstart,ieqend,i4d,cn,ci,cc,cplus,where,match);
    }
    Loop(iter,0,iterproc) 
    { 
      Loop(i,0,noproc)
      {   
        if (nameproc[i]=='c') eqncenter(ieqstart,ieqend,cplus,rhs,ci, ch);
        else if (nameproc[i]==ctdma[0]) tdmasolve(tdma1[0],istdma[0],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[1]) tdmasolve(tdma1[1],istdma[1],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[2]) tdmasolve(tdma1[2],istdma[2],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[3]) tdmasolve(tdma1[3],istdma[3],ci,cplus,rhs,ch);
        else { printout("normal"," procedure not found\n"); if (noproc==1) break; }
        
        eqnrhs(ieqstart,ieqend,cn,ci,cc,ch,rhsin,  rhs); 
      }
      eqnerrorwho(ieqstart,ieqend,rhs,cplus,ny, err,&imin,&imax); /* error analysis */ 
      if (iter==0 || iter==iterproc-1 || err[0]>err[3])
      { 
        iexpand(whois[imin],i4d,iimin); iexpand(whois[imax],i4d,iimax);
        printout("normal"," iter %d  err-rms %lg    err-min %lg at %d %d %d %d    err-max %lg at %d %d %d %d\n",iter,
                err[0],err[1],iimin[0],iimin[1],iimin[2],iimin[3],
                err[2],iimax[0],iimax[1],iimax[2],iimax[3]);			
      }
      err[3]=err[0];
    }
  }
  /* end of solve , clean up*/
  printout("normal","end solve\n"); 
  if (camin>0)   /* reinstate center coef */
    Loop(i,ieqstart,ieqend)
    if (cn[i]>0) cc[i][0]=csave[i];
  if (nomit>0) /* reinstart + sign of coef_n  */
  {Loop(i,ieqstart,ieqend) if (cn[i]<0) cn[i]=-cn[i]; }
  Loop (k,0,4) if (tdma1[k]>0) tdmafree(tdma1[k]);
  
  /* set matched points */
  Loop(i,0,npts[0]) ch[i]=ch[match[i]];
}

