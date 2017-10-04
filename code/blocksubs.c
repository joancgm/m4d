/* contains c_blocksetabc, blockcoef, blockrhs, blocktdma, blocktop  */
#include "global.h"
#include "iexpand.h"
#include "tdmasubs.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void coefprint1(int i, int *cn, int **ci, double **cc, int jrep, FILE *fprint);

/* multigrid solver for pressure equations, divide up p-eqn ijk so far no 4th dimension */

/* set blockijk? to be used by blocksubs used by eqnsolvep for muntigrid solve */
/*  blockijk= ibtot[3] number of blocks in each direction 
 bl[3][ibtot[i]+1] start ip location for each block, with next value giving end of loop
 ipinb[3][ip4d[i]] the ib location for each ip point
 set up so far to omit boundary ip=0 and ip=last equations */	
void c_blocksetabc(FILE *fpin, FILE *fprint)
{  
  int *i4d; double *a[3]; char *csym; /* need */
  char label;int in[3]; double *ain; /* read */
  int *blockijk; /*create */
  int *ibtot,*bl[3],*ipinb[3];
  int i,j,je,k;
  char name[]="blockijkx";
  
  i4d=(int *)need("idim4d");
  csym=(char *)need("csym");
  a[0]=(double *)need("abcd"); a[1]=a[0]+i4d[0]; a[2]=a[1]+i4d[1];
  
  ain=(double *)tmalloca(max(i4d[0],max(i4d[1],i4d[2])),'d');
  
  label=read1charname(fpin);
  Loop(i,0,3) { in[i]=readint(fpin); in[i]=min(in[i],i4d[i]); }
  
  name[8]=label;
  printout("normal"," creating %s for %d by %d by %d blocks\n",name,in[0]-1,in[1]-1,in[2]-1);
  /* check input including check for 1-d repeat which may fail */
  Loop(i,0,3)
  {
    if (in[i]<2)
    {
      printout("error c_blocksetabc","error must have at least 2 points each direction %d\n",i);
      exitm4d(0);
    }  
    if (in[i]==2 && in[(i+1)%3]==2 && csym[6+2*(i+2)%3]=='r')
    {
      printout("error c_blocksetabc","1-d solve repeat not allowed - it can blow up\n");
      exitm4d(0);
    }
  }
  blockijk=(int *)createarray(name,3+in[0]+in[1]+in[2]
                              +(i4d[0]+1)+(i4d[1]+1)+(i4d[2]+1),'i',0); /* ib,bl,ipinb, */
  
  Loop(i,0,3) blockijk[i]=in[i]-1;
  ibtot=blockijk;
  bl[0]=blockijk+3; bl[1]=bl[0]+in[0]; bl[2]=bl[1]+in[1];
  ipinb[0]=bl[2]+in[2]; ipinb[1]=ipinb[0]+i4d[0]+1;  ipinb[2]=ipinb[1]+i4d[1]+1; 
  
  Loop(i,0,3) 
  {  
    if (in[i]==i4d[i])  /* use individual lines do not read a b or c's */
    {  
      printout("normal"," individual lines for  %d \n",i);
      Loop(k,0,in[i]) bl[i][k]=k+1;
      ipinb[i][0]=-1; Loop(k,1,i4d[i]) ipinb[i][k]=k-1; ipinb[i][i4d[i]]=-1;
    }
    else
    {  
      printout("normal"," read a b or c %d\n",in[i]);
      Loop(k,0,in[i])     /* read a b or c's */
      { 
        ain[k]=readdouble(fpin); 
        printout("normal"," %lg",ain[k]); 
        if (k>0 && ain[k]<ain[k-1]) 
        { 
          printout("error c_blocksetabc"," input error \n"); 
          exitm4d(0); 
        }
      }
      ain[0]=max(ain[0],a[i][0]);
      ain[in[i]-1]=min(ain[in[i]-1],a[i][i4d[i]-1]);
      printout("normal","\n");
      Loop(k,0,i4d[i]+1) ipinb[i][k]=-1;
      j=0; je=1;
      Loop(k,0,i4d[i])  
      {  
        if (a[i][k]<ain[0]) continue;
        if (a[i][k]==ain[0]) bl[i][0]=k+1;
        if (a[i][k]==ain[in[i]-1]) bl[i][je]=k+1;
        if (a[i][k]>=ain[in[i]-1]) break;
        if (a[i][k]>=ain[je]) {	j++; je=min(in[i]-1,je+1); bl[i][j]=k+1;}
        ipinb[i][k+1]=j;
      }
    }
    printout("normal","bl%d",i); Loop(k,0,in[i]) printout("normal"," %d",bl[i][k]); printout("normal","\n");
    printout("normal","ipinb%d",i); Loop(k,0,i4d[i]+1) printout("normal"," %d",ipinb[i][k]); printout("normal","\n");
    
  }
}

/* returned : acnb, acib, aneqsb; accb, acplus */
int blockcoef(char m, int *cn, int **ci, double **cc, int *i4dp, int *wherep, char *cltp, int **acnb, int ***acib, double ***accb, int *aneqsb, double **acplus)
{  
  int *cnb,**cib; double **ccb,*cplus;
  int i,j,neqsb,ipt,ieqb,ieq,iptbcen,iptb,ipta;
  int ipp[4],ipb[4]={0,0,0,1},ip[4]={0,0,0,0},ib[4]={0,0,0,0};
  
  int *cibb,ny,nyomit,ni;
  double *ccbb; 
  
  int *blockijk; 
  int *ibtot,*bl[3],*ipinb[3]; 
  char name[]="blockijkx";
  int ieqbug=-1; 
  
  name[8]=m;
  blockijk=(int *)find(name);
  if (blockijk==0) return 0;
  printout("normal","blockcoef %c\n",m);
  ibtot=blockijk; 
  bl[0]=blockijk+3; bl[1]=bl[0]+ibtot[0]+1; bl[2]=bl[1]+ibtot[1]+1;
  ipinb[0]=bl[2]+ibtot[2]+1; ipinb[1]=ipinb[0]+i4dp[0];  ipinb[2]=ipinb[1]+i4dp[1]; 
  
  Loop(i,0,3) ipb[i]=ibtot[i];
  printout("normal"," ipb %d %d %d %d\n",ipb[0],ipb[1],ipb[2],ipb[3]);
  neqsb=ipb[0]*ipb[1]*ipb[2];
  
  cnb=(int *)tmalloca(neqsb,'i');
  Loop(i,0,neqsb) cnb[i]=0;
  cib=(int **)tmalloca(neqsb,'p');
  ccb=(double **)tmalloca(neqsb,'p');
  cplus=(double *)tmalloca(neqsb,'d');
  
  cibb=(int *)tmalloca(neqsb,'i');
  ccbb=(double *)tmalloca(neqsb,'d');
  
  Loop(ib[0],0,ipb[0]) Loop(ib[1],0,ipb[1]) Loop(ib[2],0,ipb[2]) /* each cont group */
  {    
    ieqb=ib[0]+ipb[0]*(ib[1]+ipb[1]*ib[2]);
    nyomit=0;
    ni=0; 
    iptbcen=In4(ib,ipb);
    if (ieqb==ieqbug) printout("normal","ib0 %d ip0 %d %d ib1 %d ip1 %d %d ib2 %d ip2 %d %d iptbc %d\n",ib[0],bl[0][ib[0]],bl[0][ib[0]+1],ib[1],bl[1][ib[1]],bl[1][ib[1]+1],ib[2],bl[2][ib[2]],bl[2][ib[2]+1],iptbcen);
    Loop(ip[0],bl[0][ib[0]],bl[0][ib[0]+1])
    Loop(ip[1],bl[1][ib[1]],bl[1][ib[1]+1])
    Loop(ip[2],bl[2][ib[2]],bl[2][ib[2]+1])    /* each equation in group */
    {
      ipta=In4(ip,i4dp);
      if (cltp[ipta]=='o') nyomit=1;
      ieq=wherep[ipta];
      if (ieqb==ieqbug) printout("normal","ip %d %d %d ieq %d \n",ip[0],ip[1],ip[2],ieq);
      if (ieq<0) continue;
      if (ieqb==ieqbug) coefprint1(ieq,cn,ci,cc,1,stdout); 
      if (cn[ieq]==0) continue;
      
      Loop(i,0,cn[ieq])
      { 
        ipt=ci[ieq][i];
        iexpand(ipt,i4dp,ipp);
        if (ipinb[0][ipp[0]]<0) continue;
        if (ipinb[1][ipp[1]]<0) continue;
        if (ipinb[2][ipp[2]]<0) continue;   /* omit for points which are not considered in blocks */
        iptb=ipinb[0][ipp[0]]+ibtot[0]*(ipinb[1][ipp[1]]+ibtot[1]*ipinb[2][ipp[2]]);
        if (ieqb==ieqbug) printout("normal","ipp %d %d %d ipinb %d %d %d iptb %d cc %lg\n",
                                 ipp[0],ipp[1],ipp[2],ipinb[0][ipp[0]],ipinb[1][ipp[1]],ipinb[2][ipp[2]]
                                 ,iptb,cc[ieq][i]);
        if (ni==0)
        {        /* set as +- including repeat bndry values */
          Loop(j,0,7) { cibb[j]=iptbcen; ccbb[j]=0; }
          if (ib[0]>0) cibb[1]--; else cibb[1]+=ibtot[0]-1;
          if (ib[0]<ibtot[0]-1) cibb[2]++; else cibb[2]-=ibtot[0]-1;
          if (ib[1]>0) cibb[3]-=ipb[0]; else cibb[3]+=ipb[0]*(ibtot[1]-1);
          if (ib[1]<ibtot[1]-1) cibb[4]+=ipb[0]; else cibb[4]-=ipb[0]*(ibtot[1]-1);
          if (ib[2]>0) cibb[5]-=ipb[0]*ipb[1]; else cibb[5]+=ipb[0]*ipb[1]*(ibtot[2]-1);
          if (ib[2]<ibtot[2]-1) cibb[6]+=ipb[0]*ipb[1]; else cibb[6]-=ipb[0]*ipb[1]*(ibtot[2]-1);
          ni=7;
          if (ieqb==ieqbug) {printout("normal","cibb"); Loop(j,0,7) printout("normal"," %d",cibb[j]); printout("normal","\n"); }
        }
        ny=0;
        Loop(j,0,ni) 
        if (iptb==cibb[j]) 
        { 
          ccbb[j]+=cc[ieq][i]; 
          ny=1;
          if (ieqb==ieqbug)  printout("normal","j %d ccbb %lg\n",j,ccbb[j]);
          break; 
        }
        if (ny==0) 
        { 
          cibb[ni]=iptb; 
          ccbb[ni]=cc[ieq][i]; 
          if (ieqb==ieqbug)  printout("normal","ni %d ccbb %lg\n",ni,ccbb[ni]);
          ni++; 
        }
      }
    }  /* ip */
    cnb[ieqb]=ni;
    cplus[ieqb]=0;
    
    if (ni>0)
    {  
      cib[ieqb]=(int *)tmalloca(ni,'i');
      ccb[ieqb]=(double *)tmalloca(ni,'d');
      Loop(i,0,ni) { cib[ieqb][i]=cibb[i]; ccb[ieqb][i]=ccbb[i]; }
      
      Loop(i,0,cnb[ieqb]) if (ccb[ieqb][i]>0) cplus[ieqb]+= ccb[ieqb][i]; 
    }
    if (ieqb==ieqbug)
    {  
      printout("normal","ib %d %d %d cplus %lg\n",ib[0],ib[1],ib[2],cplus[ieqb]);
      coefprint1(ieqb,cnb,cib,ccb,1,stdout); 
    }
    if (nyomit==1 && cnb[ieqb]>0) cnb[ieqb]=-cnb[ieqb];
  }	
  acnb[0]=cnb; acib[0]=cib; accb[0]=ccb; aneqsb[0]=neqsb; acplus[0]=cplus;
  printout("normal","  %d equations\n",neqsb); 
  return 1;
}
void blockrhs(char m, int *i4dp, int *wherep, double *rhs, int *cnb, double **arhsb)
{    
  double *rhsb;
  char name[]="blockijkx";
  int *blockijk,*ibtot,*bl[3];
  int ib[3],ip[4]={0,0,0,0},ipta,ieqb,ieq;
  double rhsbmin,rhsbmax,rhsbrms; int neqb,i;
  int ieqbug=-1;
  
  name[8]=m;
  blockijk=(int *)need(name); 
  
  /* printout("normal","blockrhs\n"); */
  ibtot=blockijk; 
  bl[0]=blockijk+3; bl[1]=bl[0]+ibtot[0]+1; bl[2]=bl[1]+ibtot[1]+1;
  rhsb=arhsb[0];
  if (rhsb==0)
  {  
    rhsb=(double *)tmalloca(ibtot[0]*ibtot[1]*ibtot[2],'d');
    arhsb[0]=rhsb;
  }
  Loop(ib[2],0,ibtot[2]) Loop(ib[1],0,ibtot[1]) Loop(ib[0],0,ibtot[0])  /* each cont group */
  {
    ieqb=ib[0]+ibtot[0]*(ib[1]+ibtot[1]*ib[2]);
    rhsb[ieqb]=0;
    
    Loop(ip[0],bl[0][ib[0]],bl[0][ib[0]+1])
    Loop(ip[1],bl[1][ib[1]],bl[1][ib[1]+1])
    Loop(ip[2],bl[2][ib[2]],bl[2][ib[2]+1])    /* each equation in group */
    {
      ipta=In4(ip,i4dp);
      ieq=wherep[ipta];
      if (ieq<0) continue;
      rhsb[ieqb]+=rhs[ieq];
      if (ieqb==ieqbug) printout("normal","ip %d %d %d rhs %lg rhsb %lg\n",ip[0],ip[1],ip[2],rhs[ieq],rhsb[ieqb]);
    }
    /* if (cnb[ieqb]<=0) rhsb[ieqb]=0; */
    /* if (ib[0]==0) printout("normal","\n %d %d",ib[1],ib[2]);
     printout("normal"," %lg",rhsb[ieqb]);
     */
  }
  neqb=ibtot[0]*ibtot[1]*ibtot[2];
  rhsbmin=rhsb[0]; rhsbmax=rhsb[0]; 
  rhsbrms=rhsb[0]*rhsb[0];
  Loop(i,1,neqb)
  { 
    rhsbmin=min(rhsbmin,rhsb[i]); rhsbmax=max(rhsbmax,rhsb[i]);
    rhsbrms+=rhsb[i]*rhsb[i];
  }
  rhsbrms=sqrt(rhsbrms/(double)neqb);
  printout("normal"," rhsb  rhs %lg min %lg max %lg\n",rhsbrms,rhsbmin,rhsbmax); 
}
/* two pass ijk solve, update rhs and dp after each call */
void blocktdma(char m, double *cplus, int *cnb, int **cib, double **ccb, double *rhsb, double **adpb)
{
  int i,j,n,n1,n2,maxd,ib[3],ieq,iptb,neqb;
  double *dx,*cpn,*bbn,*ccn,*cmn,d,*rhs,*dpb;
  double dpbmin,dpbmax,rhsbmin,rhsbmax,rhsbrms; 
  int iter;
  char name[]="blockijkx"; int *ibtot;
  int bug=0; 
  /* int ieqbug=-1;  */
  
  name[8]=m;
  ibtot=(int *)need(name);
  maxd=max(ibtot[0],max(ibtot[1],ibtot[2]));
  neqb=ibtot[0]*ibtot[1]*ibtot[2];
  
  if (bug > 0) printout("normal","blocktdma maxd %d, neqb %d \n",maxd,neqb);	
  /* Loop(i,0,neqb) { printout("normal"," %d",cnb[i]); if (i%100==0) printout("normal","\n"); } 
  printout("normal","\n"); */
  dpb=adpb[0]; 
  if (dpb==0)
  {  
    dpb=(double *)tmalloca(neqb,'d');
    adpb[0]=dpb;
  }
  Loop(i,0,neqb) dpb[i]=0;
  
  dx=(double *)smalloca(maxd*5,'d');
  cpn=dx+maxd; bbn=cpn+maxd; ccn=bbn+maxd; cmn=ccn+maxd;
  rhs=(double *)smalloca(neqb,'d');
  Loop(i,0,neqb) if (cnb[i]>0) rhs[i]=rhsb[i]; else rhs[i]=0;
  
  if (bug>0) printout("normal","space allocated\n");
  Loop(iter,0,2)
  Loop(n,0,3)
  {  
    n1=(n+1)%3; n2=(n+2)%3;
    if (bug > 0) printout("normal","tdma n n1 n2 %d %d %d iter %d\n",n,n1,n2,iter);  
    Loop(ib[n1],0,ibtot[n1]) Loop(ib[n2],0,ibtot[n2])
    {  
      if (bug>1 && n==2) printout("normal"," ib1 %d ib2 %d\n",ib[n1],ib[n2]); 
      Loop(ib[n],0,ibtot[n])
      { 
        ieq=ib[0]+ibtot[0]*(ib[1]+ibtot[1]*ib[2]);
        if (bug>1 && n==2) printout("normal","ieq %d cnb %d cplus %lg\n",ieq,cnb[ieq],cplus[ieq]);
        if (cnb[ieq]<=0 || cplus[ieq]<=0) 
        { 
          ccn[ib[n]]=1; cmn[ib[n]]=0; cpn[ib[n]]=0; bbn[ib[n]]=0;
        }
        else
        { ccn[ib[n]]=cplus[ieq];
          cmn[ib[n]]=ccb[ieq][1+2*n];
          cpn[ib[n]]=ccb[ieq][2+2*n];
          bbn[ib[n]]=rhs[ieq];
        }
        if (bug>1 && n==2)	printout("normal"," %d cc %lg cm %lg cp %lg bb %lg\n",
                                  ib[n],ccn[ib[n]],cmn[ib[n]],cpn[ib[n]],bbn[ib[n]]); 
      }
      
      if (ibtot[n]==1) dx[0]=bbn[0]/ccn[0];
      else if (cmn[0]!=0 || cpn[ibtot[n]-1]!=0) /* repeat tdma */
        tdmar(cmn,ccn,cpn,bbn,dx,ibtot[n]);
      else /* reg tdma solve */
      {
        if (ibtot[n]==2)
        {
          dx[1] = (bbn[1]-cmn[1]*bbn[0]/ccn[0])/(ccn[1]-cmn[1]*cpn[0]/ccn[0]);
          dx[0] = bbn[0]/ccn[0] -dx[1]*cpn[0]/ccn[0];
        }
        else   /* max>2 */
        {
          cpn[0] = -cpn[0]/ccn[0];
          bbn[0] = bbn[0]/ccn[0];
          for (j=1;j<ibtot[n]-1;j++)
          { 
            d = ccn[j]+cmn[j]*cpn[j-1];
            bbn[j] = (bbn[j]-cmn[j]*bbn[j-1])/d;
            cpn[j] = -cpn[j]/d;
          }
          dx[ibtot[n]-1] = (bbn[ibtot[n]-1]-cmn[ibtot[n]-1]*bbn[ibtot[n]-2])
          /(ccn[ibtot[n]-1]+cmn[ibtot[n]-1]*cpn[ibtot[n]-2]);
          for (j=ibtot[n]-2;j>=0;j--) dx[j] = bbn[j]+cpn[j]*dx[j+1];
        }
      }
      Loop(ib[n],0,ibtot[n])
      { 
        iptb=ib[0]+ibtot[0]*(ib[1]+ibtot[1]*ib[2]);
        dpb[iptb]+=dx[ib[n]];
        if (bug>1 && n==2)	printout("normal"," %lg",dx[ib[n]]); 
      }
      if (bug>1 && n==2) printout("normal","\n"); 
    }
    Loop(i,0,neqb)
    {
      if (cnb[i]>0)
      {  
        rhs[i]=rhsb[i];
        Loop(j,0,cnb[i]) rhs[i]-=ccb[i][j]*dpb[cib[i][j]];
      }
    }
    rhsbmin=rhs[0]; rhsbmax=rhs[0]; dpbmin=dpb[0]; dpbmax=dpb[0];
    rhsbrms=rhs[0]*rhs[0];
    Loop(i,1,neqb)
    { 
      rhsbmin=min(rhsbmin,rhs[i]); rhsbmax=max(rhsbmax,rhs[i]);
      dpbmin=min(dpbmin,dpb[i]); dpbmax=max(dpbmax,dpb[i]);
      rhsbrms+=rhs[i]*rhs[i];
    }
    rhsbrms=sqrt(rhsbrms/(double)neqb);
    if (bug > 0) printout("normal"," rhsb rhs %lg %lg %lg dpb %lg %lg\n",
           rhsbrms,rhsbmin,rhsbmax,dpbmin,dpbmax); 
  }
  /* Loop(ib[2],0,ibtot[2]) Loop(ib[1],0,ibtot[1]) Loop(ib[0],0,ibtot[0])
	{  iptb=ib[0]+ibtot[0]*(ib[1]+ibtot[1]*ib[2]);
   if (ib[0]==0) printout("normal","\n %d %d",ib[1],ib[2]);
   printout("normal"," %lg",rhs[iptb]);
   if (iptb==ieqbug) printout("normal","*");
	}
	*/
  if (bug > 0) printout("normal"," blocktdma done\n");	 
  free(dx); free(rhs);
  
}
/* ------------------------- */
void blocktop(char m, int *i4dp, double *dpb, int *cnb, int **cib, double *dp)
{   
  int *blockijk; 
  int *ibtot,*bl[3],*ipinb[3],iptb; 
  char name[]="blockijkx";
  int ip[4]={0,0,0,0};
  /*int i,ieqbug=-1; */
  
  name[8]=m;
  blockijk=(int *)need(name);
  
  /* printout("normal","blocktop\n"); */
  ibtot=blockijk; 
  bl[0]=blockijk+3; bl[1]=bl[0]+ibtot[0]+1; bl[2]=bl[1]+ibtot[1]+1;
  
  ipinb[0]=bl[2]+ibtot[2]+1; ipinb[1]=ipinb[0]+i4dp[0];  ipinb[2]=ipinb[1]+i4dp[1]; 
  
  Loop(ip[2],0,i4dp[2]) Loop(ip[1],0,i4dp[1]) Loop(ip[0],0,i4dp[0])
  {  
    if (ipinb[0][ip[0]]<0) continue;
    if (ipinb[1][ip[1]]<0) continue;
    if (ipinb[2][ip[2]]<0) continue;   /* omit for points which are not considered in blocks */
    iptb=ipinb[0][ip[0]]+ibtot[0]*(ipinb[1][ip[1]]+ibtot[1]*ipinb[2][ip[2]]);
    dp[In4(ip,i4dp)] += dpb[iptb];
    /* if (ieqbug>=0 && cnb[ieqbug]>0)
     Loop(i,0,cnb[ieqbug]) 
     if (iptb==cib[ieqbug][i]) printout("normal","ip %d %d %d %d ib %d %d %d  %d dpb %lg dp %lg \n",
     ip[0],ip[1],ip[2],In4(ip,i4dp),
     ipinb[0][ip[0]],ipinb[1][ip[1]],ipinb[2][ip[2]]
     ,iptb,dpb[iptb],	dp[In4(ip,i4dp)]);
     */
  }
}

