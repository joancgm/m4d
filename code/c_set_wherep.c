/* contains c_set_wherep */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* set number of independent p pts and wherep in point list  (-1 if not) */
/* also set cltp    - indices indicating c.v. types:
 f-flow F-nearwall flow side, s-solid, S-near wall solid side
 (flow is anything other than clt=w or s) */

void c_set_wherep(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*matchpc; char *clt;   /* needed arrays */
  int *noindppts,*wherep,*whoisp;  char *cltp; /* set arrays */
  
  int i,j,iallp,nop,i4dp[4],idiv,it,ip[4],ii[4],ia[3],iw,is;
  char cc;
  int two[3]={2,2,2};
  int bug=0;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  idiv=iallp/i4d[3];
  matchpc=(int *)need("matchpc");
  clt=(char *)need("clt");
  /*   cltp */
  cltp=(char *)createarray("cltp",iallp,'c',0);
  Loop(i,0,iallp)
  { 
    iexpand(i,i4dp,ip);
	 iw=0; is=0; ii[3]=ip[3];
	 Loop3(ia,0,two)
    { 
      Loop(j,0,3) ii[j]=min(i4d[j]-1,max(ip[j]-1+ia[j],0));
      cc=clt[In4(ii,i4d)];
      if (cc=='w') iw++;
      else if (cc=='s') is++;
    }
	 if (iw==0 && is==0) cltp[i]='f';
	 else if (iw>0 && is==0) cltp[i]='F';
	 else if (iw>0 && is>0) cltp[i]='S';
	 else if (iw==0 && is>0) cltp[i]='s';
  }
  if (bug>0) printout("normal","cltp set\n");
  /*---------- noindppts, wherep, whoisp ------------  */
  noindppts=(int *)createarray("noindppts",i4d[3]+2,'i',0); 
  Loop(i,0,i4d[3]+2) noindppts[i]=0;
  wherep=(int *)createarray("wherep",iallp,'i',0);
  Loop(i,0,iallp) wherep[i]=-1;  /* initialize as matched point */
  /* printout("normal","t1\n"); */
  nop=0;
  Loop(i,0,iallp) 	  
  {
    /* printout("normal"," %d %d %d",i,matchpc[i]); */
    if (matchpc[i]<0) 
    { 
      wherep[i]=nop; 
      nop++;  
      it=i/idiv;
      noindppts[it+2]=nop;
    }
  }
  noindppts[0]=nop;
  printout("normal","nop %d\n",nop);
  whoisp=(int *)createarray("whoisp",noindppts[0],'i',0);
  nop=0;
  Loop(i,0,iallp) if (matchpc[i]<0) {whoisp[nop]=i; nop++; }
  
  printout("normal","set wherep,  whoisp, noindppts = ");
  Loop(i,0,i4d[3]+2) printout("normal","  %d",noindppts[i]);
  printout("normal","\n");
}
