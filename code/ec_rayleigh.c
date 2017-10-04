/* contains c_rayleigh */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* 
 
 */

void ec_rayleigh(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep;   /* need */
  double *a[4],*u[3],*w,*walldist; 
  char *cltp;
  double *ray, *elpd, *vtpd;  /* set */
  double ase[3][2]; int iline;  /* read */
  double cyw,cyi,cypd,cely;
  
  int i,j,k,i4dm[4],iallm,iall,ip[4],i4dp[4],ii[4],mmid,mm,mmidm,mmidp,ipe;
  int iis[4],iie[4],id1,id2;
  int jwmax=0,jumin,jumax,jc7[2],jvisc[2];
  double f,wmax,umin,umax,els=0,vtpda=0,umag,ywmax,elsp;
  
  i4d=(int *)need("idim4d"); 
  a[0]=(double *) need("abcd");
  Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  ii[3]=0; ip[3]=0,iis[3]=0;
  iallm=Prod4(i4dm);
  iall=Prod4(i4d);
  
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  w=(double *)need("wrate");
  walldist=(double *)need("walldist");
  
  /* input */
  Loop(i,0,3) 
  {
    Loop(j,0,2) ase[i][j]=readdouble(fpin);
    iis[i]=findex(ase[i][0],a[i],i4d[i],&f);
    if (f>.9) iis[i]++;
    iie[i]=findex(ase[i][1],a[i],i4d[i],&f);
    if (f>.9) iie[i]++;
    printout("normal","n %d, an=%lg, %lg, i[n]= %d %d\n",
            i,ase[i][0],ase[i][1],iis[i],iie[i]);
  }
  iline=readint(fpin);
  cyw=readdouble(fpin);
  cyi=readdouble(fpin);
  cypd=readdouble(fpin);
  cely=readdouble(fpin);
  printout("normal","\niline= %d, cwy %lg, cyi %lg, cypd %lg, cely %lg\n",
          iline,cyw,cyi,cypd,cely);
  iline--;
  id1=(iline+1)%3;
  id2=(iline+2)%3;
  printout("normal","il 1 2 3 %d %d %d\n",iline,id1,id2);
  
  ray=(double *)createarray("ray",iallm,'d',0);
  vtpd=(double *)createarray("vtpd",iallm,'d',0);
  elpd=(double *)createarray("elpd",iallm,'d',0);
  Loop(i,0,iallm) { ray[i]=0; vtpd[i]=0; elpd[i]=0; }
  
  /* loop over midlines */
/*  printout("normal"," n+1  n+2  jw yw W jmin Umin jmax Umax Ls vtpd jc7 jvtpd\n"); */
  Loop(ii[id1],iis[id1],iie[id1])  Loop(ii[id2],iis[id2],iie[id2]) 
  {
    /* find umin, umax */
    ii[iline]=iis[iline];
    mmid=In4(ii,i4dm);
    mm=In4(ii,i4d);
  /*  printout("normal","ijk %d %d %d mmid mm %d %d\n",ii[0],ii[1],ii[2],mmid,mm); */
    jumin=ii[iline]; 
    umag=sqrt(u[0][mm]*u[0][mm]+u[1][mm]*u[1][mm]+u[2][mm]*u[2][mm]);
    umin=umag; 
    jumax=jumin; 
    umax=umin;
  /*  printout("normal","umag %lg w %lg j's %d %d %d\n",umax,w[mmid],jwmax,jumin,jumax); */
    
    Loop(ii[iline],iis[iline],iie[iline])
    {
      mm=In4(ii,i4d);
  /*    printout("normal","ijk %d %d %d mmid mm %d %d",ii[0],ii[1],ii[2],mmid,mm); */
      umag=sqrt(u[0][mm]*u[0][mm]+u[1][mm]*u[1][mm]+u[2][mm]*u[2][mm]);
      if (umag<umin) { umin=umag; jumin=ii[iline]; }
      if (umag>umax) { umax=umag; jumax=ii[iline]; }
  /*    printout("normal","umag %lg w %lg j's %d %d %d\n",umax,w[mmid],jwmax,jumin,jumax); */
    }
    printout("normal"," %d:%d %d:%d j Umax  %d %lg \n",id1,ii[id1],id2,ii[id2],jumax,umax);
    /* find wmax locations */
    Loop(k,iis[iline]+1,iie[iline]-2)
    {
      ii[iline]=k;
      ip[0]=ii[0]+1; ip[1]=ii[1]+1; ip[2]=ii[2]+1;
      ipe=wherep[In4(ip,i4dp)];
      if (ipe<0) continue;
      mmid=In4(ii,i4dm);
      ii[iline]--; mmidm=In4(ii,i4dm);
      ii[iline]+=2; mmidp=In4(ii,i4dm);
      ii[iline]--;
      if (w[mmid]==0) continue;      
     /* printout("normal"," j 3w's %d %lg %lg %lg \n",k,w[mmidm],w[mmid],w[mmidp]); */
      if (w[mmid]>=w[mmidm] && w[mmid]>=w[mmidp] && 2*w[mmid]>w[mmidm]+w[mmidp])
      {
        els=(umax-umin)/w[mmid];
      /*  printout("normal"," y,crit %lg %lg\n",ywallm[mmid],cyw*els); */
        if (walldist[ipe]<cyw*els) continue;
        /* found wmax location */
        vtpda=els*els*w[mmid];
        jwmax=k;
        ywmax=walldist[ipe];
        wmax=w[mmid];
        printout("normal"," %d:%d %d:%d j wmax elsv vpda %d %lg %lg %lg",
               id1,ii[id1],id2,ii[id2],jwmax,wmax,els,vtpda);
        jc7[0]=iie[iline]-1; jc7[1]=iis[iline]; jvisc[0]=iie[iline]-1; jvisc[1]=iis[iline];
        Loop(ii[iline],iis[iline],iie[iline]-1)
        { 
          mmid=In4(ii,i4dm);
          ip[0]=ii[0]+1; ip[1]=ii[1]+1; ip[2]=ii[2]+1;
          ipe=wherep[In4(ip,i4dp)];
          if (ipe<0) continue;
          if(abs(ywmax-walldist[ipe]) < cyi*els)
          { 
            if (ii[iline]<jc7[0]) jc7[0]=ii[iline];
            if (ii[iline]>jc7[1]) jc7[1]=ii[iline];
            ray[mmid]=1;
          }
          if( abs(ywmax-walldist[ipe]) < cypd*els)
          { 
            if (ii[iline]<jvisc[0]) jvisc[0]=ii[iline];
            if (ii[iline]>jvisc[1]) jvisc[1]=ii[iline];
            {
              elsp=els;
              if (cely>0) elsp=min(els,cely*walldist[ipe]);
              if (elpd[mmid]==0 || elpd[mmid]>elsp) 
              {
                elpd[mmid]=elsp;
                vtpd[mmid]=elpd[mmid]*elpd[mmid]*wmax;
              }
            }
          }
        }
        printout("normal"," jL %d %d jv %d %d\n",jc7[0],jc7[1],jvisc[0],jvisc[1]);
      }
    }
  }
}