/* contains c_interp_ptog, interp_ptog */
#include "global.h"

/* interpolate from p points to grid points  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Value3d(xv,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*xv[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*xv[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*xv[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*xv[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*xv[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*xv[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*xv[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*xv[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])

void interp_ptog(int nprop, double **pp, double **pg);

void c_interp_ptog(FILE *fpin, FILE *fprint)
{ 
  int *i4d; /* need */
  double *pp[20], *pg[20];  /* from and to */
  int n; char*name,*named; /* input parmaters */
  int i,j,iall;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  
  n=readint(fpin);  if (n>20) n=20;
  printout("normal","interp");
  Loop(i,0,n) 
  {
    name=readname(fpin);
    pp[i]=(double *)need(name);
    named=readname(fpin);
    printout("normal",", %s to %s",name,named);
    pg[i]=(double *)find(named);
    if (pg[i]==0) 
    { 
      pg[i]=(double*)createarray(named,iall,'d',0);
      Loop(j,0,iall) pg[i][j]=0;
    }
  }
  printout("normal","\n");
  interp_ptog(n,pp,pg);
}
/* -------------------------- */
void  interp_ptog(int nprop, double **pp, double **pg)
{ 
  int *i4d,**whoelse,*match; double *xyz,*xpc;  /* need */
  int *itrange;  /* use if available */
  
  int i,k,n,m,i4dp[4],iallp,iall,its,ite;
  int id[4],ilo[4],ihi[4],ipt,iappt;
  double f[3],*x[3],xx[3],*ppi;
  double tol=.00001;
  int iter;
  
  int itermax=11;
  double fmin,fmax,delf,df[3],xt[3],d[11],dx[3],dxdk[3][3],fa[3],fb[3],det,xtt[3],dold;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  iall=Prod4(i4d);
  whoelse=(int **)need("whoelse");
  match=(int *)need("match");
  xyz=(double *)need("xyz");
  xpc=(double *)need("xyzp");
  itrange=(int *)find("itrange");
  if (itrange>0) { its=itrange[0]; ite=itrange[1]; }
  else {its=0; ite=i4d[3]-1;}
  
  Loop(id[3],its,ite+1) 
  {  
    iappt=id[3]*i4dp[0]*i4dp[1]*i4dp[2]; 
    x[0]=xpc+iappt; 
    x[1]=x[0]+iallp; x[2]=x[1]+iallp;
    /*  printout("normal","id[3] %d iappt %d iapdt %d \n",id[3],iappt,iapdt); */
    Loop3(id,0,i4d)
    {
      ipt=In4(id,i4d);
      if (match[ipt] != ipt) continue;
      Loop(i,0,3) xx[i]=xyz[ipt+iall*i];
      
      /* new version with old as backup */      {
        /* printout("normal","id %d %d %d\n",id[0],id[1],id[2]); */
        fmin=.1*tol;
        fmax=1-.1*tol;
        delf=fmax-fmin;
        Loop(i,0,3) 
        {
          ilo[i]=id[i]; 
          ihi[i]=id[i]+1;
          f[i]=.5;
          if (ilo[i]==0) f[i]=fmin;
          else if (ilo[i]==i4d[i]-1) f[i]=fmax;
        }
        Loop(i,0,3)
        {
          xt[i]=Value3d(x[i],i4dp,ilo,f);
          dx[i]=xx[i]-xt[i];
        }
        d[0]=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
        Loop(iter,1,itermax)
        {
          Loop(k,0,3)
          { 
            Loop(n,0,3) {fa[n]=f[n],fb[n]=f[n]; }
            fa[k]=fmin;
            fb[k]=fmax;
            Loop(i,0,3) dxdk[i][k]=Value3d(x[i],i4dp,ilo,fb)-Value3d(x[i],i4dp,ilo,fa);
          }
          det=dxdk[0][0]*(dxdk[1][1]*dxdk[2][2]-dxdk[1][2]*dxdk[2][1])
          +dxdk[0][1]*(dxdk[1][2]*dxdk[2][0]-dxdk[1][0]*dxdk[2][2])
          +dxdk[0][2]*(dxdk[1][0]*dxdk[2][1]-dxdk[1][1]*dxdk[2][0]);
          /* printout("normal","iter %d f %lg %lg %lg det %lg\n",iter,f[0],f[1],f[2],det); */
          if (det==0) break;
          Loop(i,0,3) dx[i]=xx[i]-xt[i];
          df[0]=dx[0]*(dxdk[1][1]*dxdk[2][2]-dxdk[1][2]*dxdk[2][1])
          +dxdk[0][1]*(dxdk[1][2]*dx[2]-dx[1]*dxdk[2][2])
          +dxdk[0][2]*(dx[1]*dxdk[2][1]-dxdk[1][1]*dx[2]);
          df[1]=dxdk[0][0]*(dx[1]*dxdk[2][2]-dxdk[1][2]*dx[2])
          +dx[0]*(dxdk[1][2]*dxdk[2][0]-dxdk[1][0]*dxdk[2][2])
          +dxdk[0][2]*(dxdk[1][0]*dx[2]-dx[1]*dxdk[2][0]);
          df[2]=dxdk[0][0]*(dxdk[1][1]*dx[2]-dx[1]*dxdk[2][1])
          +dxdk[0][1]*(dx[1]*dxdk[2][0]-dxdk[1][0]*dx[2])
          +dx[0]*(dxdk[1][0]*dxdk[2][1]-dxdk[1][1]*dxdk[2][0]);
          Loop(i,0,3) 
          {
            df[i]/=(det*delf);
            fa[i]=max(fmin,min(fmax,f[i]+df[i]));
            if (ilo[i]==0) fa[i]=fmin;
            else if (ilo[i]==i4d[i]-1) fa[i]=fmax;
            df[i]=fa[i]-f[i];
          }
          /* printout("normal","df %lg %lg %lg\n",df[0],df[1],df[2]); */
          if (abs(df[0])<tol && abs(df[1])<tol && abs(df[2])<tol) break;
          
          Loop(k,0,5)
          { 
            Loop(i,0,3) fa[i]=f[i]+df[i];
            Loop(i,0,3)
            {
              xtt[i]=Value3d(x[i],i4dp,ilo,fa);
              dx[i]=xx[i]-xtt[i];
            }
            d[iter]=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]; 
            /* printout("normal","d %lg dm1 %lg\n",d[iter],d[iter-1]); */
            if (d[iter]==0 || d[iter]<d[iter-1])
            { 
              Loop(i,0,3) {f[i]=fa[i]; xt[i]=xtt[i]; }
              break;
            }
            if (k<4) Loop(i,0,3) df[i]*=.2;
            d[iter]=d[iter-1];
          } /* k loop */
          /* printout("normal","k=%d\n",k); */
          if (d[iter]==d[iter-1]) /* use old method */
          {
            /* printout("normal","new %d %d %d f %lg %lg %lg d d0 %lg %lg\n",
             id[0],id[1],id[2],f[0],f[1],f[2],d[iter],d[0]); */
            Loop(i,0,3) {fa[i]=f[i]; f[i]+=ilo[i]; }
            k=findex3d(xx,x,i4dp,ilo,ihi,tol,f);
            Loop(i,0,3) f[i]-=ilo[i];
            Loop(i,0,3)
            {
              xtt[i]=Value3d(x[i],i4dp,ilo,f);
              dx[i]=xx[i]-xtt[i];
            }
            Loop(i,0,3) df[i]=f[i]-fa[i];
            dold=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
            /*  printout("normal","old method k %d f= %lg %lg %lg dold %lg \n",k,f[0],f[1],f[2],dold); */
            if (dold>d[iter]) Loop(i,0,3) f[i]=fa[i];
            break;
          }
        } /* iter */
        
        Loop(i,0,3) 
        {
          if (ilo[i]==0) f[i]=0;
          else if (ilo[i]==i4d[i]-1) f[i]=1;
        }
        /* printout("normal","finterp %lg %lg %lg\n",f[0],f[1],f[2]); */
        Loop(i,0,nprop)
        {
          ppi=pp[i]+iappt;
          pg[i][ipt]=Value3d(ppi,i4dp,ilo,f);
          if  (whoelse[ipt]>0)
          { 
            k=whoelse[ipt][0];
            Loop(m,1,k+1) pg[i][whoelse[ipt][m]]=pg[i][ipt];
          }
        }
      } 
    }
  }
}

