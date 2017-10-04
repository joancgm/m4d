/* contains c_interp_ptod, interp_ptod */
#include "global.h"

/* interpolate from p points to double points */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void interp_ptod(int nprop, double **pp, double **pd);

void c_interp_ptod(FILE *fpin, FILE *fprint)
{ 
  int *i4d; /* need */
  double *pp[20], *pd[20];  /* from and to */
  int n; char*name,*named; /* input parmaters */
  int i,j,nod;
  
  i4d=(int *)need("idim4d");
  nod=(2*i4d[0]-1)*(2*i4d[1]-1)*(2*i4d[2]-1)*i4d[3];
  
  n=readint(fpin);  if (n>20) n=20;
  printout("normal","interp");
  Loop(i,0,n) 
  {
    name=readname(fpin);
    pp[i]=(double *)need(name);
    named=readname(fpin);
    printout("normal",", %s to %s",name,named);
    pd[i]=(double *)find(named);
    if (pd[i]==0) 
    { 
      pd[i]=(double*)createarray(named,nod,'d',0);
      Loop(j,0,nod) pd[i][j]=0;
    }
  }
  printout("normal","\n");
  
  interp_ptod(n,pp,pd);
}
/* ---------------------------- */
void  interp_ptod(int nprop, double **pp, double **pd)
{ 
  int *i4d; double *xd,*xpc;  /* need */
  int *itrange;   /* use if available */
  
  int i,i4dp[4],idd[4],iallp,ialld,its,ite;
  int id[4],ilo[3],ihi[3],ipt,iappt,iapdt;
  double f[3],*x[3],xx[3];
  double tol=.001;
  int iter;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  Loop(i,0,3) idd[i]=2*i4d[i]-1; idd[3]=i4d[3];
  ialld=Prod4(idd);
  xd=(double *)need("xyzdouble");
  xpc=(double *)need("xyzp");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1]; }
  else {its=0; ite=i4d[3]-1;}
  
  Loop(id[3],its,ite+1) 
  { 
    iappt=id[3]*i4dp[0]*i4dp[1]*i4dp[2]; 
    iapdt=id[3]*idd[0]*idd[1]*idd[2]; 
    x[0]=xpc+iappt; 
    x[1]=x[0]+iallp; x[2]=x[1]+iallp;
    /*  printout("normal","id[3] %d iappt %d iapdt %d \n",id[3],iappt,iapdt); */
    Loop3(id,0,idd)
    {
      Loop(i,0,3)
      { if (id[i]==0) {ilo[i]=0; ihi[i]=0; }
      else if (id[i]==idd[i]-1) {ilo[i]=id[i]/2+1; ihi[i]=ilo[i]; }
      else 
      {
        ilo[i]=(id[i]+1)/2; 
        if (id[i]%2==1) ihi[i]=ilo[i];  
        else ihi[i]=ilo[i]+1;
      }
      }
      /*  printout("normal","id %d %d %d ilo %d %d %d ihi %d %d %d",
       id[0],id[1],id[2],ilo[0],ilo[1],ilo[2],ihi[0],ihi[1],ihi[2]); */
      ipt=In4(id,idd);
      Loop(i,0,3) xx[i]=xd[ipt+ialld*i];
      Loop(i,0,3) { f[i]=ilo[i]+ihi[i]; f[i]*=.5; }
      /*  printout("normal","  ipt %d xx %lg %lg %lg  f %lg %lg %lg \n",ipt,xx[0],xx[1],xx[2],f[0],f[1],f[2]); */
      iter=findex3d(xx,x,i4dp,ilo,ihi,tol,f);
      /*   printout("normal","iter %d f %lg %lg %lg\n",iter,f[0],f[1],f[2]); */
      Loop(i,0,nprop)
      { 
        interp3d(pp[i]+iappt,i4dp,f,pd[i]+iapdt+ipt);
      }
    }
  }
}

