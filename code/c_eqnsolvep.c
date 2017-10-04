/* contains c_eqnsolvep */
#include "global.h"
#include "tdmasubs.h"
#include "eqnsubs.h"
#include "blocksubs.h"
#include "psolve1d.h"
#include "psleepcalc.h"
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* solve pressure point  equations 
 read information on type of solve to be used 
 any point which lack coefficients (or has all zero ones) will be left unchanged  */
void c_eqnsolvep(FILE *fpin, FILE *fprint)
{
  int *i4d,*npts,*neqs,*where,*whois,*cn,**ci,*match; /* needed arrays, see below for real names  */
  double **cc,*rhsin,*cpsleep;  char *cltp;
  int *itrange,*ipfix,*ipeqexit, *cnall; /* use if available */
  double *ch,*eqnerr;  /* calculate */
  
  char *nameproc; int iterproc;   /* read procedure name and iterations */
  
  int i,k=0,iter,ns,ne,i4dp[4],noproc,imin,iimin[4],imax,iimax[4];
  double *rhs,*cplus,errinit[3],err[4]={0,0,0,0};
  Tdma *tdma1[4]={0,0,0,0};
  int istdma[4]={0,0,0,0}; char ctdma[]="tijk";
  char *ny=0;
  int iblock=0,nyset,iblockmax=5,**cib[5],neqsb[5],*cnb[5]={0,0,0,0,0}; 
  char nameblock[]="_____"; 
  double *cplusb[5],**ccb[5],*rhsb[5]={0,0,0,0,0},*dpb[5]={0,0,0,0,0};
  double pfix;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  npts=(int *)need("noppts");
  neqs=(int*)need("noindppts");
  where=(int *)need("wherep");
  whois=(int *)need("whoisp");
  match=(int *)need("matchpc");
  cnall=(int *)need("cpc_n"); /* for bndry eqns and/or no cpc1_n */
  cn=(int *)find("cpc1_n");
  if (cn==0) cn=cnall;
  ci=(int **)need("cpc_i");
  cc=(double **)need("cpc_c");
  rhsin=(double *)need("rhsc");
  cltp=(char *)need("cltp");
  itrange=(int *)find("itrange");
  ipfix=(int *)find("ipfix");
  ipeqexit=(int *)find("ipeqexit");
  if (itrange==0) {ns=0; ne=neqs[0]; }
  else {ns=neqs[itrange[0]+1]; ne=neqs[itrange[1]+2]; }
  istdma[3]=ns;
  
  ch=(double *)createarray("dp",npts[0],'d',0);
  eqnerr=(double *)createarray("eqnerrp",6,'d',0);
  Loop(i,0,npts[0]) ch[i]=0;
  
  /* initial calcs and analysis */
  rhs=(double *)tmalloca(neqs[0],'d');
  cplus=(double *)tmalloca(neqs[0],'d');
  /* omit eqns with cltp='o' by changing sigh of c, */
  Loop(i,ns,ne) if (cltp[whois[i]]=='o' && cn[i]>0) cn[i]=-cn[i];
  /* also for eqns listed in ipeqexit */
  if (ipeqexit>0)
    Loop(i,0,ipeqexit[0]) if (cn[ipeqexit[i+1]]>0) cn[ipeqexit[i+1]]=-cn[ipeqexit[i+1]];
  eqncplus(ns,ne,cn,cc, cplus); /* cplus=sum of + coef */
  /*  printout("normal","cplus set \n"); */
  eqnrhs(ns,ne,cnall,ci,cc,ch,rhsin,  rhs);  /* set rhs with  ch included */
  /*  printout("normal"," rhs set \n"); */
  eqnerrorwho(ns,ne,rhs,cplus,ny, err,&imin,&imax); /* error analysis */ 
  iexpand(whois[imin],i4dp,iimin); iexpand(whois[imax],i4dp,iimax);
  printout("normal"," init:  perr-rms %lg    perr-min %lg at %d %d %d %d    perr-max %lg at %d %d %d %d\n",
          err[0],err[1],iimin[0],iimin[1],iimin[2],iimin[3],
          err[2],iimax[0],iimax[1],iimax[2],iimax[3]);
  Loop(i,0,3) errinit[i]=err[i];
  eqnerr[0]=err[0];
  eqnerr[1]=max(-err[1],err[2]);
  printout("normal","enter solvep\n");
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
      if (nameproc[i]=='b' || nameproc[i]=='c' || nameproc[i]=='I' ||
          nameproc[i]=='J' || nameproc[i]=='K') continue;
      if (nameproc[i]=='t') 
      { if (tdma1[0]==0) tdma1[0]=tdmafromcoef(ns,ne,cn,ci,cc,cplus,where); }
      else if (nameproc[i]=='i') 
      { if (tdma1[1]==0) tdma1[1]=tdmaijkcoef(0,ns,ne,i4dp,cn,ci,cc,cplus,where,match); }
      else if (nameproc[i]=='j') 
      {  if(tdma1[2]==0) tdma1[2]=tdmaijkcoef(1,ns,ne,i4dp,cn,ci,cc,cplus,where,match); }
      else if (nameproc[i]=='k')
      {  if (tdma1[3]==0)  tdma1[3]=tdmaijkcoef(2,ns,ne,i4dp,cn,ci,cc,cplus,where,match); }
      else  /* abc blocks */
      {  
        nyset=0;
        Loop(k,0,iblock) if (nameproc[i]==nameblock[k]) nyset=1;
        if (nyset==1) continue;
        if (iblock==iblockmax)
        { 
          printout("error c_eqnsolvep"," redimesion for more blocks %d\n",iblockmax); 
          exitm4d(0); 
        }
        k=iblock;
        nyset=blockcoef(nameproc[i],cn,ci,cc,i4dp,where,cltp,&cnb[k],&cib[k],&ccb[k],
                        &neqsb[k],&cplusb[k]);
        if (nyset>0) { iblock++; nameblock[k]=nameproc[i];	}
        else 
        { 
          printout("wrror c_eqnsolvep"," info for block %c not found\n",nameproc[i]); 
          exitm4d(0); 
        }
      }
      
    }
    Loop(iter,0,iterproc) 
    {
      Loop(i,0,noproc)
      {  
        k=-1;
        if (nameproc[i]=='c')  eqncenter(ns,ne,cplus,rhs,ci, ch);
        else if (nameproc[i]=='b')  eqnbndry(ipeqexit,cnall,ci,cc,rhsin,rhs,ch);
        else if (nameproc[i]==ctdma[0]) tdmasolve(tdma1[0],istdma[0],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[1]) tdmasolve(tdma1[1],istdma[1],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[2]) tdmasolve(tdma1[2],istdma[2],ci,cplus,rhs,ch);
        else if (nameproc[i]==ctdma[3]) tdmasolve(tdma1[3],istdma[3],ci,cplus,rhs,ch);
        else if (nameproc[i]=='I') psolve1d(0,cn,ci,cc,i4dp,where,cltp,ipfix,rhs,ch);
        else if (nameproc[i]=='J') psolve1d(1,cn,ci,cc,i4dp,where,cltp,ipfix,rhs,ch);
        else if (nameproc[i]=='K') psolve1d(2,cn,ci,cc,i4dp,where,cltp,ipfix,rhs,ch);
        else 
          Loop(k,0,iblock) if (nameproc[i]==nameblock[k]) 
          {  
            blockrhs(nameproc[i],i4dp,where,rhs,cnb[k], &rhsb[k]);
            blocktdma(nameproc[i],cplusb[k],cnb[k],cib[k],ccb[k],rhsb[k],&dpb[k]);
            blocktop(nameproc[i],i4dp,dpb[k],cnb[k],cib[k], ch);
            break;
          }
        
        if (k==iblock) 
        { 
          printout("normal"," procedure %c not found\n",nameproc[i]); 
          if (noproc==1) break; 
        }
        
        eqnrhs(ns,ne,cnall,ci,cc,ch,rhsin,  rhs); 
      }
      eqnerrorwho(ns,ne,rhs,cplus,ny, err,&imin,&imax); 
      if (iter==0 || iter==iterproc-1 || err[0]>err[3])
      { 
        iexpand(whois[imin],i4dp,iimin); iexpand(whois[imax],i4dp,iimax);
        printout("normal"," iter %d  perr-rms %lg    perr-min %lg    perr-max %lg\n",
                iter,	 err[0],err[1],err[2]);
      }
      err[3]=err[0];
    }
  }
  eqnerr[2]=err[0];
  eqnerr[3]=max(-err[1],err[2]);
  if (eqnerr[0]>0) eqnerr[4]=eqnerr[2]/eqnerr[0]; else eqnerr[4]=0;
  if (eqnerr[1]>0) eqnerr[5]=eqnerr[3]/eqnerr[1]; else eqnerr[5]=0;
  /* fill in changes for sleeeping points */
  cpsleep=(double *)need("cpsleep"); 
  psleepcalc(fprint,ch,cpsleep);
  /* clean up */
  Loop(i,ns,ne) if (cltp[whois[i]]=='o' && cn[i]<0) cn[i]=-cn[i];
  if (ipeqexit>0)
    Loop(i,0,ipeqexit[0]) if (cn[ipeqexit[i+1]]<0) cn[ipeqexit[i+1]]=-cn[ipeqexit[i+1]];

  if (ipfix>0) /* set dp 0 at this point and relatively change all others */
  { 
    pfix=ch[In4(ipfix,i4dp)];
    Loop(i,0,npts[0]) ch[i]-=pfix;
  }
  printout("normal","end solve\n"); 
  Loop (k,0,4) if (tdma1[k]>0) tdmafree(tdma1[k]);
}
