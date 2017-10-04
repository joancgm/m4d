/* contains c_inletinit, c_contcpcinlet, c_inletreset */
#include "global.h"
/* subs for flow in/out boundary conditions */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* set up for multiple  inflow regions with total pressure/velocity ratio  fixed */
/* if all three velocity ratios are zero, set velocity ratio from current solution */
/* regions with fixed velocity, independent of pressure do not need to be added to cont eqns */

void c_inletinit(FILE *fpin, FILE *fprint)
{ 
  int *i4d; char *clt;/* need */
  int *flowinn; double *flowinprop;  /* create */
  
  int i,j,k,itot,itott,is[4],ie[4],jst,nin,L,ip,ii[4],nyset=0,ipt,idin[4];
  double *propin; /* temporary */
  double *pg=0,*U[3],*rho=0; /* needed if setting UL/Upot from current solution */
  double ptmax,*pt,*urat[3],upot;
  /* flowinprop: ptmax, U1,U2,U3/|Upotential|   */
  
  i4d=(int *)need("idim4d"); 
  clt=(char *)need("clt");
  nin=readint(fpin); 
  printout("normal"," %d inlet regions\n",nin);
  /* flowinn: 0 number of regions, then for each region, itot,ist,iend,L,dipL */
  flowinn=(int *)createarray("flowinn",1+5*nin,'i',0); 
  propin=(double *)tmalloca(4*nin,'d');
  flowinn[0]=nin;    /* read number of regions */
  itott=0;
  Loop(i,0,nin)       /* read 4-d indices for regions */
  {
    itot=1;
    Loop(j,0,4)  
    { 
      is[j]=readint(fpin)-1;   /* later tag to abc  ?? */ 
      is[j]=max(is[j],0);
      ie[j]=readint(fpin)-1; 
      ie[j]=min(ie[j],i4d[j]-1);
      itot*=(ie[j]+1-is[j]);
    }
    printout("normal","is,ie=%d %d,   %d %d,   %d  %d,   %d  %d",
            is[0],ie[0],is[1],ie[1],is[2],ie[2],is[3],ie[3]);
    flowinn[1+i*5]=itot;    
    flowinn[2+i*5]=In4(is,i4d);
    flowinn[3+i*5]=In4(ie,i4d);
    L=-1;
    Loop(j,0,3) if (ie[j]==is[j]) {L=j; break;}
    printout("normal","   L=%d  ",L);
    if (L==-1) 
    {
      printout("error c_inletinit","Error inlet surface direction\n");
      exitm4d(0); 
    }
    flowinn[4+i*5]=L;
    if (is[L]==0) ip=0;
    else if (is[L]==i4d[L]-1) ip=-1;
    else /* check clt in midrange */
    {
      ip=0;
      Loop(j,0,4) ii[j]=(is[j]+ie[j])/2;  
      ii[L]=is[L]+1;
      if (clt[In4(ii,i4d)]=='w' ||  clt[In4(ii,i4d)]=='s') ip=-1;
    }
    flowinn[5+i*5]=ip;
    if (ip==0) printout("normal"," p-side +\n");
    else if (ip==-1) printout("normal"," p-side -\n");
    itott+=itot;
    
    /* read constants, separate for each region: ptmax, U1,U2,U3/|Upotential */
    Loop(j,0,4) propin[j+4*i]=readdouble(fpin);
    if (propin[1+4*i]==0 && propin[2+4*i]==0 && propin[3+4*i]==0)
    {
      printout("normal","ptmax, UL/Upotential set from currect solution for region\n");
      nyset=1;
    }
    else printout("normal","ptmax= %lg, UL/Upotential= %lg, %lg, %lg for region\n",
                 propin[0+4*i],propin[1+4*i],propin[2+4*i],propin[3+4*i]);
  }
  /* printout("normal","%d inflow boundary points\n",itott); */
  flowinprop=(double *)createarray("flowinprop",itott*4,'d',0);
  if (nyset==1) 	 /* needed to set from current */
  {
    rho=(double *)need("rho");
    pg=(double *)need("pg");
    U[0]=(double *)need("U1");
    U[1]=(double *)need("U2");
    U[2]=(double *)need("U3");
  }
  jst=0;
  Loop(i,0,nin)
  { 
    Loop(k,0,flowinn[1+5*i])
    Loop(j,0,4) flowinprop[jst+k+j*flowinn[1+5*i]]=propin[j+4*i];
    if (propin[1+4*i]==0 && propin[2+4*i]==0 && propin[3+4*i]==0)
    { 
      /* set from current */
      iexpand(flowinn[2+i*5],i4d,is);
      iexpand(flowinn[3+i*5],i4d,ie); Loop(j,0,4) ie[j]++;
      ipt=In4(is,i4d);
      ptmax=pg[ipt]+.5*rho[ipt]*
      (U[0][ipt]*U[0][ipt]+U[1][ipt]*U[1][ipt]+U[2][ipt]*U[2][ipt]);
      Loop(ii[3],is[3],ie[3]) 
      Loop(ii[2],is[2],ie[2]) Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      {
        ipt=In4(ii,i4d);
        ptmax=max(ptmax,pg[ipt]+.5*rho[ipt]*
                  (U[0][ipt]*U[0][ipt]+U[1][ipt]*U[1][ipt]+U[2][ipt]*U[2][ipt]));
      } 
      Loop(j,0,4) idin[j]=ie[j]-is[j];
      pt=flowinprop+jst; 
      Loop(j,0,3) urat[j]=pt+(j+1)*flowinn[1+5*i];
      Loop(ii[3],is[3],ie[3]) 
      Loop(ii[2],is[2],ie[2]) Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      { 
        k=ii[0]-is[0]+idin[0]*(ii[1]-is[1]+idin[1]*(ii[2]-is[2]+idin[2]*(ii[3]-is[3])));
        ipt=In4(ii,i4d);
        upot=sqrt((ptmax-pg[ipt])/(.5*rho[ipt]));
        pt[k]=ptmax;
        Loop(L,0,3) urat[L][k]=U[L][ipt]/upot;
      }
    }
    jst+=4*flowinn[1+5*i];
  }
}
/*-----------------------------------------*/
void c_contcpcinlet(FILE *fpin, FILE *fprint)
{
  int *i4d, *wherep;  double *pp,*rho;  /* need */
  int *flowinn; double *flowinprop;  /* need */
  int *itrange; /* use if available */
  int *cpc_n, **cpc_i;  double **cpc_c;  /* modify */
  
  int i,j,k,jst,is[4],ie[4],ii[4],ipt,idin[4],ic[4],ia,ib,i4dp[4],ip[4],L,L2,L3,ieq;
  double *pt,*urat[3],area[4][3],fda,pta,rhom,upot;
  int bug=0;
  
  flowinn=(int *)find("flowinn");
  if (flowinn==0) return;
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  itrange=(int *)find("itrange");
  wherep=(int *)need("wherep"); 
  flowinprop=(double *)need("flowinprop");
  cpc_n=(int*)need("cpc_n"); 
  cpc_i=(int**)need("cpc_i"); 
  cpc_c=(double**)need("cpc_c"); 
  rho=(double *)need("rho");
  pp=(double *)need("pp");
  geomcinit();
  
  printout("normal"," %d inflow regions:\n",flowinn[0]);
  jst=0;
  Loop(i,0,flowinn[0])
  { 
    iexpand(flowinn[2+i*5],i4d,is);
    iexpand(flowinn[3+i*5],i4d,ie);  Loop(j,0,4) ie[j]++;
    printout("normal","i=%d <%d, j=%d <%d, k=%d <%d, t=%d <%d\n",
            is[0],ie[0],is[1],ie[1],is[2],ie[2],is[3],ie[3]);
    Loop(j,0,4) idin[j]=ie[j]-is[j];
    L=flowinn[4+i*5];
    Loop(j,0,3) if (j!=L) ie[j]--;
    L2=(L+1)%3; L3=(L+2)%3;
    pt=flowinprop+jst; 
    Loop(j,0,3) urat[j]=pt+(j+1)*flowinn[1+5*i];
    Loop(ii[3],is[3],ie[3]) 
    {
      if (itrange>0) { if (ii[3]<itrange[0]|| ii[3]>itrange[1]) continue; }
      Loop(ii[2],is[2],ie[2]) Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      {
        geomcavector(area,ii,L);
        if (bug>0)	  printout("normal","ii %d %d %d %d L %d\n",ii[0],ii[1],ii[2],ii[3],L);
        fda=0; pta=0; rhom=0;
        Loop(ia,0,2) Loop(ib,0,2)
        {
          Loop(j,0,4) ic[j]=ii[j];
          ic[L2]+=ia; ic[L3]+=ib;
          ipt=In4(ic,i4d);
          k=ic[0]-is[0]+idin[0]*(ic[1]-is[1]+idin[1]*(ic[2]-is[2]+idin[2]*(ic[3]-is[3])));
          if (bug>0)	  printout("normal","ia %d ib %d ic %d %d %d %d k %d\n",ia,ib,ic[0],ic[1],ic[2],ic[3],k);
          fda+=area[ia+2*ib][0]*urat[0][k]+area[ia+2*ib][1]*urat[1][k]
          +area[ia+2*ib][2]*urat[2][k];
          pta+=.25*pt[k];   /* approximate */
          rhom+=.25*rho[ipt];
        }
        Loop(j,0,3) ip[j]=ii[j]+1; ip[3]=ii[3];
        ip[L]+=flowinn[5+i*5];
        if (bug>0)	  printout("normal","fda %lg pt %lg rhom %lg ip %d %d %d %d\n",
                               fda,pta,rhom,ip[0],ip[1],ip[2],ip[3]);
        if ((fda>0 && flowinn[5+i*5]!=0) || (fda<0 && flowinn[5+i*5]!=-1))
        {
          printout("normal","bad inlet flow direction, omit for  ip= %d %d %d %d\n",
                  ip[0],ip[1],ip[2],ip[3]); 
          continue;
        }
        ipt=In4(ip,i4dp);
        ieq=wherep[ipt];
        if (ieq<0) continue;
        if (cpc_n[ieq]==0) continue;
        if (cpc_i[ieq][0]!=ipt) 
        {
          printout("normal","bad pc eqn, inlet omitted for ip= %d %d %d %d\n",
                  ip[0],ip[1],ip[2],ip[3]); 
          continue;
        }
        if (pta<pp[ipt])
        {
          printout("normal","pp %lg > pt %lg, inlet omitted for ip= %d %d %d %d\n",
                  pp[ipt],pta, ip[0],ip[1],ip[2],ip[3]); 
          continue;
        }
        upot=sqrt(2.*(pta-pp[ipt])/rhom);  /* incompressible flow */
        if (bug>0) 
          printout("normal","upot %lg cpc %lg dcpc %lg\n",upot,cpc_c[ieq][0],fda/upot);
        cpc_c[ieq][0]+=fda/upot;
      }
    }
    /*	  printout("normal", "  pt=%lg, U1/Up= %lg, U2/Up=%lg, U3/Up=%lg\n",   
     flowinprop[jst+0*flowinn[1+i*5]], flowinprop[jst+1*flowinn[1+i*5]], 
     flowinprop[jst+2*flowinn[1+i*5]],flowinprop[jst+3*flowinn[1+i*5]]); */
    jst+=4*flowinn[1+i*5];
  }
}
/* ----------------------- */
void c_inletreset(FILE *fpin, FILE *fprint)
{ 
  int *i4d;  double *rho,*xyz,*xpc;  /* need */
  int *flowinn; double *flowinprop;  /* need */
  int *itrange;  /* use if available */
  double *u[3],*pp; /* modify  u at inlet and pp to match about inlet */
  
  int i,j,k,L,ii[4],ip[4],i4dp[4],is[4],ie[4],jst,idin[4],iappt,ilo[3],ihi[3],iter,iallp,iall,ipt,iptout;
  double *pt,*urat[3],*x[3],xx[3],f[3],pg,upot,uratsq,upotoldsq;
  double tol=.001;
  
  flowinn=(int *)find("flowinn");
  if (flowinn==0) 
  {
    printout("normal","warming: omitting  c_inletreset, flowinn not available\n");
    return;
  }
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  iall=Prod4(i4d);
  xyz=(double *)need("xyz");
  xpc=(double *)need("xyzp");
  itrange=(int *)find("itrange");
  flowinprop=(double *)need("flowinprop");
  rho=(double *)need("rho");
  pp=(double *)need("pp");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  
  /* reset outside pp */
  Loop(i,0,flowinn[0])
  { 
    iexpand(flowinn[2+i*5],i4d,is);
    iexpand(flowinn[3+i*5],i4d,ie); 
    L=flowinn[4+i*5];
    ie[3]++; ie[L]++;
    Loop(j,0,3) if (j!=L)  ie[j]=min(i4d[j]+1,ie[j]+2); 
    Loop(ii[3],is[3],ie[3]) 
    {
      if (itrange>0) { if (ii[3]<itrange[0]|| ii[3]>itrange[1]) continue; }
      Loop(ii[2],is[2],ie[2]) Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      {
        Loop(j,0,3) ip[j]=ii[j]; ip[3]=ii[3];
        ip[L]+=1+flowinn[5+i*5];
        ipt=In4(ip,i4dp);
        /*	  printout("normal","set pp %lg from pt %d %d %d %d  to pt",pp[ipt],ip[0],ip[1],ip[2],ip[3]); */
        ip[L]=ii[L]+1;
        if (flowinn[5+i*5]==0) ip[L]--;
        iptout=In4(ip,i4dp);
        /*	  printout("normal"," %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]); */
        pp[iptout]=pp[ipt];
      }
    }
  }
  /* reset velocity */
  jst=0;
  Loop(i,0,flowinn[0])
  { 
    iexpand(flowinn[2+i*5],i4d,is);
    iexpand(flowinn[3+i*5],i4d,ie); Loop(j,0,4) ie[j]++;
    printout("normal","i=%d <%d, j=%d <%d, k=%d <%d, t=%d <%d\n",
            is[0],ie[0],is[1],ie[1],is[2],ie[2],is[3],ie[3]); 
    Loop(j,0,4) idin[j]=ie[j]-is[j];
    pt=flowinprop+jst; 
    Loop(j,0,3) urat[j]=pt+(j+1)*flowinn[1+5*i];
    Loop(ii[3],is[3],ie[3]) 
    {
      if (itrange>0) { if (ii[3]<itrange[0]|| ii[3]>itrange[1]) continue; }
      iappt=ii[3]*i4dp[0]*i4dp[1]*i4dp[2]; 
      x[0]=xpc+iappt; 
      x[1]=x[0]+iallp; x[2]=x[1]+iallp;
      Loop(ii[2],is[2],ie[2]) Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      { 
        k=ii[0]-is[0]+idin[0]*(ii[1]-is[1]+idin[1]*(ii[2]-is[2]+idin[2]*(ii[3]-is[3])));
        uratsq=urat[0][k]*urat[0][k]+urat[1][k]*urat[1][k]+urat[2][k];
        ipt=In4(ii,i4d);
        if (uratsq==0) 
        { 
          Loop(j,0,3) u[j][ipt]=0;
          continue;
        }
        Loop(j,0,3)	
        {
          if (ii[j]==0) {ilo[j]=0; ihi[j]=0; }
          else if (ii[j]==i4d[j]-1) {ilo[j]=ii[j]+1; ihi[j]=ilo[j]; }
          else {ilo[j]=ii[j];  ihi[j]=ii[j]+1;   }
        }
        upotoldsq=(u[0][ipt]*u[0][ipt]+u[1][ipt]*u[1][ipt]+u[2][ipt]*u[2][ipt])/uratsq;
        Loop(j,0,3) xx[j]=xyz[ipt+iall*j];
        Loop(j,0,3) { f[j]=ilo[j]+ihi[j]; f[j]*=.5; }
        iter=findex3d(xx,x,i4dp,ilo,ihi,tol,f);
        interp3d(pp+iappt,i4dp,f,&pg);
        if (pt[k]-pg<.25*upotoldsq*.5*rho[ipt])  
        {
          printout("normal","warning high p %lg,  U at %d %d %d %d halved\n"
                  ,pg,ii[0],ii[1],ii[2],ii[3]);
          Loop(j,0,3) u[j][ipt]=.5*u[j][ipt];
          continue;
        }
        upot=sqrt(2*(pt[k]-pg)/rho[ipt]);   /* incompressible flow */
        Loop(j,0,3) u[j][ipt]=upot*urat[j][k];
        /*   printout("normal","set U at %d %d %d %d to %lg %lg %lg\n",ii[0],ii[1],ii[2],ii[3],
         urat[0][k],urat[1][k],urat[2][k]); */
      }
    }
    jst+=4*flowinn[1+5*i];
  }
}


