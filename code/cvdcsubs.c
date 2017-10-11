/* contains cvdclimit1, cvdcactive, cvdccorner, cvdclines, cvdcface, cvdcfill, cvdcsleep, c_cvdcparm, c_cvdcinit, c_cvdcreset  */
#include "global.h"
#include "cmatchfollow.h"

/* control volume boundaries based  on cont c.v. center movement  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define In4low(ii,idim) (ii[0]-1+idim[0]*(ii[1]-1+idim[1]*(ii[2]-1+idim[2]* ii[3])))
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])

/* for cvdcdouble */
#define Mi(i,j,k,n) ((j)+2*i4d[1]*((k)+2*i4d[2]*((i)+i4d[0]*n)))
#define Mj(i,j,k,n) ((k)+2*i4d[2]*((i)+2*i4d[0]*((j)+i4d[1]*n)))
#define Mk(i,j,k,n) ((i)+2*i4d[0]*((j)+2*i4d[1]*((k)+i4d[2]*n)))

static double flim=0;     /* shared parameters limiting   c.v. divide  movement */
static double *flimd=0;  /* point type limitations */
static char *climd=0; 
static int limd=0;  

/* --------------------- */
void cvdclimit1(double *cvdn, int ipt, int *ipr, char *clt)
{
  int i,k,n,L,ia[3],nylimd;
  int two[3]={2,2,2};
  char cta;
  nylimd=0;
  if (limd>0) /* prelim check for different point types */
  {
    cta=clt[ipt]; 
    Loop3(ia,0,two)
    { 
      if (clt[ipt+ia[0]*ipr[0]+ia[1]*ipr[1]+ia[2]*ipr[2]]!=cta) nylimd=1;
    }
  }
  Loop(L,0,3)  /* generic limits */
  { 
    if (cvdn[L]<flim) cvdn[L]=flim; 
    else if (cvdn[L]>1.-flim) cvdn[L]=1.-flim;
  }
  if (nylimd>0)   /* pt type limits */
    Loop3(ia,0,two)
  {
    k=ipt+ia[0]*ipr[0]+ia[1]*ipr[1]+ia[2]*ipr[2];
    Loop(i,0,limd)
    {
      if (clt[k]==climd[i]) /* 1 match */
        Loop(L,0,3)
		{ 
        n=k-ia[L]*ipr[L]+(1-ia[L])*ipr[L];
        if (clt[n]==climd[i+limd])  /* 2nd match check limits */
        {
          if (ia[L]==0) 
          { if (cvdn[L]<flimd[i]) cvdn[L]=flimd[i];}
          else 
          { if (cvdn[L]>1.-flimd[i]) cvdn[L]=1.-flimd[i]; }
        }
		} /* L */
    }
  }
}
/* ---------------------------- */
void cvdcactive(int its, int ite, double *visc, double *dt, int sizedt, double t2clim, FILE *fprint)
{ 
  int *i4d,*wherep;  double *u[3],*rho,*x[3]; char *clt;  /* needed arrays */
  int *ITER; double *TIME; /* print if available */
  double *cvdc,*cvdci,*cvdcj,*cvdck;  /* modify */
  double *cvdcrms; /* set */
  
  int i,j,iall,iallp,iallm,i4dp[4],ic[4],ip[4],ig[4],ipr[4],ia[3],ipt,iptc,i4dm[4];
  int two[3]={2,2,2};
  int ich;
  double sump[3][2], suma,sumv[3],sumt,chrms[3],rhou[3],cvdn[3],cvdold[3];
  double dx[3][3],area[3][3],vol,ruarea[3],rhom,viscp,rddt,coef;
  double dcvdcmax; int ipmax[4]={0,0,0,0}, dirmax=-1;
  int dbug=0;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  ipr[0]=1; Loop(i,1,4) ipr[i]=ipr[i-1]*i4d[i-1];
  wherep=(int *)need("wherep");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  cvdc=(double *)need("cvdc");
  cvdci=(double *)need("cvdcdouble");
  cvdcj=cvdci+4*iall; cvdck=cvdcj+4*iall;
  clt=(char *)need("clt");
  x[0]=(double *)need("xyz"); 
  x[1]=x[0]+iall; 
  x[2]=x[1]+iall;
  cvdcrms=(double *)createarray("cvdcrms",5,'d',0);
  
  
  /*  printout("normal","found arrays\n"); */
  
  ich=0; Loop(i,0,3) chrms[i]=0;
  dcvdcmax=0; 
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) 
  {
    if (dbug>0) printout("normal","ip %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]); 
    Loop(i,0,3) ic[i]=ip[i]-1; 
    ic[3]=ip[3];   ig[3]=ip[3];
    iptc=In4(ic,i4d);
    if (wherep[In4(ip,i4dp)]<0) /* set sleeps to -1 for later check */
    {
      Loop(i,0,3) cvdc[iptc+i*iall]=-1;
      continue; 
    }
    /* do for valid cont c.v. */
    Loop(i,0,3) cvdold[i]=cvdc[iptc+i*iall]; /* save old and reset middle with limits */
    Loop(i,0,3)  cvdn[i]=.5;  
    cvdclimit1(cvdn,iptc,ipr,clt);  
    Loop(i,0,3) cvdc[iptc+i*iall]=cvdn[i];  /* default values to be used for zero rhoU */
    /* make coef estimates based on parallel sides */
    Loop(i,0,3) rhou[i]=0;  rhom=0;
    Loop(i,0,3) Loop(j,0,3) dx[i][j]=0;
    Loop3(ia,0,two)     
    { 
      Loop(i,0,3) ig[i]=ic[i]+ia[i];
      ipt=In4(ig,i4d);
      Loop(i,0,3) rhou[i]+=.125*rho[ipt]*u[i][ipt];
      rhom+=.125*rho[ipt];
      Loop(j,0,3) Loop(i,0,3) dx[j][i]+=.25*x[i][ipt]*(2*ia[j]-1);
    }
    viscp=0; rddt=0;
    if (visc>0) viscp=visc[In4(ic,i4dm)];
    if (dt>0) rddt=rhom/dt[ip[3]];
    if (sizedt==iallm) rddt=rhom/dt[In4(ic,i4dm)];
    if (dbug>0) printout("normal","rhou %lg %lg %lg rhom %lg viscp %lg rddt %lg\n",
                       rhou[0],rhou[1],rhou[2],rhom,viscp,rddt);
    if (dbug>0) printout("normal","dxi %lg %lg %lg dxj %lg %lg %lg dxk %lg %lg %lg\n",
                       dx[0][0],dx[0][1],dx[0][2],dx[1][0],dx[1][1],dx[1][2],dx[2][0],dx[2][1],dx[2][2]);
    ipt=In4(ic,i4d);
    if (Sqr(rhou)>0)  /* non-zero rhoU determine optimum breakup of cont c.v. */
    {
      Loop(i,0,3) {Cross(dx[(i+1)%3],dx[(i+2)%3],area[i]); }
      if (dbug>0) printout("normal","area0 %lg %lg %lg area1 %lg %lg %lg area2 %lg %lg %lg\n",
                         area[0][0],area[0][1],area[0][2],area[1][0],area[1][1],
                         area[1][2],area[2][0],area[2][1],area[2][2]);
      vol=Dot(dx[0],area[0]);
      if (vol<=0)  /* grid error, print and stop */
      { 
        printout("error cvdcactive","grid error for cont c.v. ip=%d %d %d %d, vol= %lg\n",
                ip[0],ip[1],ip[2],ip[3],vol);
        printout("error","dxi %lg %lg %lg dxj %lg %lg %lg dxk %lg %lg %lg\n",
                dx[0][0],dx[0][1],dx[0][2],dx[1][0],dx[1][1],dx[1][2],dx[2][0],dx[2][1],dx[2][2]);
        exitm4d(0);
      }
      Loop(i,0,3) ruarea[i]=Dot(area[i],rhou);
      if (dbug>0) printout("normal","vol %lg ruarea %lg %lg %lg\n",vol,ruarea[0],ruarea[1],ruarea[2]);
      Loop(i,0,3) {sump[i][1]=ruarea[i]; sump[i][0]=-ruarea[i]; }
      if (dbug>0) printout("normal","convsump %lg %lg %lg %lg %lg %lg\n",sump[0][0],sump[0][1],
                         sump[1][0],sump[1][1],sump[2][0],sump[2][1]);
      suma=0;
      Loop3(ia,0,two) 
      { 
        coef=.25*(sump[0][ia[0]]+sump[1][ia[1]]+sump[2][ia[2]]);
        if (coef>0) suma+=coef;	
      }
      Loop(i,0,3) sumv[i]=2*viscp*Sqr(area[i])/vol;
      sumt=min(t2clim*suma,3*rddt*vol);
      if (dbug>0) printout("normal","convsum %lg +viscsum %lg %lg %lg +time sum %lg\n",
                         suma,sumv[0],sumv[1],sumv[2],sumt);
		
      if (suma>0) /* calc cvd based on sums */
      { 
        Loop(i,0,3)  cvdn[i]=.5*(1.+(-sump[i][1]+sump[i][0])/(suma+sumv[i]+sumt));
        if (dbug>0) printout("normal","cvdc %lg %lg %lg\n",cvdn[0],cvdn[1],cvdc[2]);
        cvdclimit1(cvdn,ipt,ipr,clt);
        ich++;
        Loop(i,0,3)  /* change analysis */
        {
          chrms[i]+=(cvdn[i]-cvdold[i])*(cvdn[i]-cvdold[i]);
          if (abs(cvdn[i]-cvdold[i])>dcvdcmax)
          {
            dcvdcmax=abs(cvdn[i]-cvdold[i]);
            dirmax=i; ipmax[0]=ip[0]; ipmax[1]=ip[1]; ipmax[2]=ip[2]; ipmax[3]=ip[3];
          }
        }
      }
    } /* end if rhou>0 */
    Loop(i,0,3) cvdc[ipt+i*iall]=cvdn[i];
    cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,2*(ip[2]-1)+1,ip[3])]=cvdn[0]; 
    cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdn[1]; 
    cvdck[Mk(2*(ip[0]-1)+1,2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdn[2]; 
  }
  
  if (ich>0)
  { 
    Loop(i,0,3) chrms[i]=sqrt(chrms[i]/ich);
    printout("normal","cvd_change %d pts, rms: %lg %lg %lg",ich,chrms[0],chrms[1],chrms[2]);
    printout("normal"," max %lg dir %d at %d %d %d %d",dcvdcmax,dirmax,ipmax[0],ipmax[1],ipmax[2],ipmax[3]);
    ITER=(int *)find("ITER");
    if (ITER!=0) printout("normal"," ITER %d  ",*ITER);
    TIME=(double *)find("TIME");
    if (TIME!=0) printout("normal"," TIME %lg ",*TIME);
    printout("normal","\n");
    Loop(i,0,3) cvdcrms[i]=chrms[i];
    cvdcrms[3]=dcvdcmax;
    cvdcrms[4]=dirmax;
  }
}
/* ----------------------------- */
/* fill in lines */
void cvdclines(int its, int ite, FILE *fprint)
{ 
  int *i4d,*match; char *csym;/* needed */
  double *cvdci,*cvdcj,*cvdck;  /* modify */
  int iall,it,ia,is,i,j,L,ig[4],igp[4],iga[4],igpa[4],La[2],igb[3];
  int im,nlist,nfaces,ncheck,*list,*faces,nykeep;
  int jpanic;
  double cvdcn,cvdca;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  match=(int *)need("match");
  csym=(char *)need("csym");
  cvdci=(double *)need("cvdcdouble");
  cvdcj=cvdci+4*iall; cvdck=cvdcj+4*iall;
  /* set up temporary arrays */
  im=max(i4d[0]*i4d[1],max(i4d[1]*i4d[2],i4d[2]*i4d[0]));
  list=(int *)tmalloca(im*2,'i');
  faces=(int *)tmalloca(im*2,'i');		
  
  Loop(it,its,ite+1) Loop(L,0,3)
  {
    ig[3]=it; igp[3]=it; iga[3]=it; igpa[3]=it;
    La[0]=(L+1)%3; La[1]=(L+2)%3;
    Loop(igb[L],0,i4d[L]-1)
    Loop(igb[La[0]],0,i4d[La[0]]) Loop(igb[La[1]],0,i4d[La[1]])
    {
      Loop(i,0,3) { ig[i]=igb[i]; igp[i]=ig[i];} igp[L]=ig[L]+1;
      /*printout("normal","igb %d %d %d L %d match %d %d\n",
       igb[0],igb[1],igb[2],L,match[In4(ig,i4d)],match[In4(igp,i4d)]); */
      if (match[In4(ig,i4d)]==match[In4(igp,i4d)])
      { 
        if (L==0) cvdci[Mi(ig[0],2*ig[1],2*ig[2],ig[3])]=.5;
        else if (L==1) cvdcj[Mj(2*ig[0],ig[1],2*ig[2],ig[3])]=.5;
        else cvdck[Mk(2*ig[0],2*ig[1],ig[2],ig[3])]=.5;
        continue;
      }
      nlist=1; list[0]=ig[La[0]]; list[1]=ig[La[1]];
      nfaces=0; ncheck=0; 
      jpanic=im;
      while (ncheck<nlist) /* check neighbors */
      { 
        jpanic--; if (jpanic<0) break; /* prevent an infinite loop if error */
        Loop(i,0,2) {ig[La[i]]=list[2*ncheck+i]; igp[La[i]]=ig[La[i]];}
        /*printout("normal","ig %d %d %d nlist %d ncheck %d match %d %d\n",
         ig[0],ig[1],ig[2],nlist,ncheck,match[In4(ig,i4d)],match[In4(igp,i4d)]); */
        
        Loop(ia,0,2) 
        { 
          Loop(i,0,3) iga[i]=ig[i];
          
          for (is=-1;is<2;is+=2)
          { 
            iga[La[ia]]=ig[La[ia]]+is; 
            /*printout("normal"," ia %d is %d iga %d %d %d",ia,is,iga[0],iga[1],iga[2]); */
            if (iga[La[ia]]<0 || iga[La[ia]]>i4d[La[ia]]-1)
            {
              if (csym[6+2*La[ia]]=='r' || csym[6+2*La[ia]]=='R')
              {
                list[2*nlist+ia]=iga[La[ia]]-is*i4d[La[ia]];
                list[2*nlist+(1-ia)]=iga[La[1-ia]];
                nykeep=1;
                Loop(i,0,nlist) 
                if (list[2*nlist]==list[2*i] && list[2*nlist]==list[2*i]) 
                { nykeep=0; break; }
                if (nykeep==0) continue;
                /*printout("normal"," rep %d %d  nlist %d ncheck %d\n",list[2*nlist],list[2*nlist+1],nlist,ncheck); */
                nlist++;
              }
              continue;
            }
            Loop(i,0,3) igpa[i]=iga[i]; igpa[L]++;
            /*printout("normal"," match %d %d %d %d\n",match[In4(ig,i4d)],match[In4(iga,i4d)]
             ,match[In4(igp,i4d)],match[In4(igpa,i4d)]); */
            if (match[In4(ig,i4d)]==match[In4(iga,i4d)] && 
                match[In4(igp,i4d)]==match[In4(igpa,i4d)])
            { 
              Loop(i,0,2) list[2*nlist+i]=iga[La[i]];
              nykeep=1;
              Loop(i,0,nlist) 
              if (list[2*nlist]==list[2*i] && list[2*nlist+1]==list[2*i+1]) { nykeep=0; break; }
              if (nykeep==1)
              {
                nlist++;
                faces[2*nfaces+ia]=2*min(iga[La[ia]],iga[La[ia]]-is)+1;
                faces[2*nfaces+(1-ia)]=2*min(iga[La[1-ia]],ig[La[1-ia]]);							
                /*printout("normal"," add iga %d %d %d nlist %d ncheck %d face %d %d\n",
                 iga[0],iga[1],iga[2],nlist,ncheck,faces[2*nfaces],faces[2*nfaces+1]); */
                nfaces++;
                
              }
            }
          }
        }
        ncheck++;
      }
      if (jpanic<0) 
      {
        printout("error c_cvdclines","inf loop in cvdclines, stop\n"); 
        exitm4d(0); 
      }
      /* determine value */
      cvdcn=0; j=0;
      /*printout("normal","cvdca"); */
      if (nfaces==0) /* real interpolation */
      { 
        Loop(ia,0,2) Loop(is,-1,1)
        {
          iga[L]=ig[L]; 
          iga[La[1-ia]]=2*ig[La[1-ia]];
          iga[La[ia]]=ig[La[ia]]+is;
          if (csym[6+2*La[ia]]=='r' || csym[6+2*La[ia]]=='R') 
            iga[La[ia]]=(iga[La[ia]]+i4d[La[ia]])%i4d[La[ia]];
          if (iga[La[ia]]<0 || iga[La[ia]]>i4d[La[ia]]-2) continue;
          iga[La[ia]]=2*iga[La[ia]]+1;
          if (L==0) cvdca=cvdci[Mi(iga[0],iga[1],iga[2],iga[3])];
          else if (L==1) cvdca=cvdcj[Mj(iga[0],iga[1],iga[2],iga[3])];
          else cvdca=cvdck[Mk(iga[0],iga[1],iga[2],iga[3])];
          /*printout("normal"," %lg",cvdca); */
          if (cvdca>0) {cvdcn+=cvdca; j++;}
        }
      }
      else Loop(i,0,nfaces)
      {  
        if (L==0) cvdca=cvdci[Mi(ig[0],faces[2*i],faces[2*i+1],ig[3])];
        else if (L==1) cvdca=cvdcj[Mj(faces[2*i+1],ig[1],faces[2*i],ig[3])];
        else cvdca=cvdck[Mk(faces[2*i],faces[2*i+1],ig[2],ig[3])];
        /* printout("normal"," %lg",cvdca); */
        if (cvdca>0) {cvdcn+=cvdca; j++;}
      }
      if (j==0) cvdcn=.5;
      else cvdcn/=(double)j;
      /*printout("normal"," j %d cvdcn %lg\n",j,cvdcn); */
      
      /* set values */
      Loop(i,0,nlist) 
      {
        if (L==0) cvdci[Mi(ig[0],2*list[2*i],2*list[2*i+1],ig[3])]=cvdcn;
        else if (L==1) cvdcj[Mj(2*list[2*i+1],ig[1],2*list[2*i],ig[3])]=cvdcn;
        else cvdck[Mk(2*list[2*i],2*list[2*i+1],ig[2],ig[3])]=cvdcn;
      }
      Loop(i,0,nfaces)
      {
        if (L==0) cvdci[Mi(ig[0],faces[2*i],faces[2*i+1],ig[3])]=cvdcn;
        else if (L==1) cvdcj[Mj(faces[2*i+1],ig[1],faces[2*i],ig[3])]=cvdcn;
        else cvdck[Mk(faces[2*i],faces[2*i+1],ig[2],ig[3])]=cvdcn;
      }
    }
  }
}
/* ------------------------------- */
/* fill in faces not yet set */
void cvdcface(int its, int ite, FILE *fprint)
{
  int *i4d; char *csym; double *cvdci,*cvdcj,*cvdck;
  int iall,i,j,k,ii,jj,kk,m,n;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  csym=(char *)need("csym");
  cvdci=(double *)need("cvdcdouble");
  cvdcj=cvdci+4*iall; cvdck=cvdcj+4*iall;
  Loop(n,its,ite+1)
  {
    Loop(j,0,i4d[1]) Loop(k,0,i4d[2]-1) Loop(i,0,i4d[0]-1)  /* i interp j face */
    {
      jj=2*j; kk=2*k+1;
      m=Mi(i,jj,kk,n);
      if (cvdci[m]>=0) continue;
      if (j==0)
      { 
        if (csym[8]=='r' || csym[8]=='R') 
        { 
          cvdci[m]=.5*(cvdci[Mi(i,1,kk,n)]+cvdci[Mi(i,2*(i4d[1]-1)-1,kk,n)]);
          cvdci[Mi(i,2*(i4d[1]-1),kk,n)]=cvdci[m];
        }
        else cvdci[m]=cvdci[Mi(i,1,kk,n)];
      }
      else if (j==i4d[1]-1) cvdci[m]=cvdci[Mi(i,2*(i4d[1]-1)-1,kk,n)];
      else cvdci[m]=.5*(cvdci[Mi(i,jj-1,kk,n)]+cvdci[Mi(i,jj+1,kk,n)]);
    }
    Loop(j,0,i4d[1]) Loop(k,0,i4d[2]-1) Loop(i,0,i4d[0]-1)  /* k interp j face */
    {
      jj=2*j; ii=2*i+1;
      m=Mk(ii,jj,k,n);
      if (cvdck[m]>=0) continue;
      if (j==0)
      { 
        if (csym[8]=='r' || csym[8]=='R') 
        { 
          cvdck[m]=.5*(cvdck[Mk(ii,1,k,n)]+cvdck[Mk(ii,2*(i4d[1]-1)-1,k,n)]);
          cvdck[Mk(ii,2*(i4d[1]-1),k,n)]=cvdck[m];
        }
        else cvdck[m]=cvdck[Mk(ii,1,k,n)];
      }
      else if (j==i4d[1]-1) cvdck[m]=cvdck[Mk(ii,2*(i4d[1]-1)-1,k,n)];
      else cvdck[m]=.5*(cvdck[Mk(ii,jj-1,k,n)]+cvdck[Mk(ii,jj+1,k,n)]);
    }
    Loop(j,0,i4d[1]-1) Loop(k,0,i4d[2]) Loop(i,0,i4d[0]-1)  /* i interp k face */
    {
      jj=2*j+1; kk=2*k;
      m=Mi(i,jj,kk,n);
      if (cvdci[m]>=0) continue;
      if (k==0)
      { 
        if (csym[10]=='r' || csym[10]=='R') 
        { 
          cvdci[m]=.5*(cvdci[Mi(i,jj,1,n)]+cvdci[Mi(i,jj,2*(i4d[2]-1)-1,n)]);
          cvdci[Mi(i,jj,2*(i4d[2]-1),n)]=cvdci[m];
        }
        else cvdci[m]=cvdci[Mi(i,jj,1,n)];
      }
      else if (k==i4d[2]-1) cvdci[m]=cvdci[Mi(i,jj,2*(i4d[2]-1)-1,n)];
      else cvdci[m]=.5*(cvdci[Mi(i,jj,kk-1,n)]+cvdci[Mi(i,jj,kk+1,n)]);
    }
    Loop(j,0,i4d[1]-1) Loop(k,0,i4d[2]) Loop(i,0,i4d[0]-1)/* j interp k face */
    {
      ii=2*i+1; kk=2*k;
      m=Mj(ii,j,kk,n);
      if (cvdcj[m]>=0) continue;
      if (k==0)
      { 
        if (csym[10]=='r' || csym[10]=='R') 
        {
          cvdcj[m]=.5*(cvdcj[Mj(ii,j,1,n)]+cvdcj[Mj(ii,j,2*(i4d[2]-1)-1,n)]);
          cvdcj[Mj(ii,j,2*(i4d[2]-1),n)]=cvdcj[m];
        }
        else cvdcj[m]=cvdcj[Mj(ii,j,1,n)];
      }
      else if (k==i4d[2]-1) cvdcj[m]=cvdcj[Mj(ii,j,2*(i4d[2]-1)-1,n)];
      else cvdcj[m]=.5*(cvdcj[Mj(ii,j,kk-1,n)]+cvdcj[Mj(ii,j,kk+1,n)]);
    }
    Loop(j,0,i4d[1]-1) Loop(k,0,i4d[2]-1) Loop(i,0,i4d[0])/* j interp i face */
    {
      ii=2*i; kk=2*k+1;
      m=Mj(ii,j,kk,n);
      if (cvdcj[m]>=0) continue;
      if (i==0)
      {
        if (csym[6]=='r' || csym[6]=='R') 
        { 
          cvdcj[m]=.5*(cvdcj[Mj(1,j,kk,n)]+cvdcj[Mj(2*(i4d[0]-1)-1,j,kk,n)]);
          cvdcj[Mj(2*(i4d[0]-1),j,kk,n)]=cvdcj[m];
        }
        else cvdcj[m]=cvdcj[Mj(1,j,kk,n)];
      }
      else if (i==i4d[0]-1) cvdcj[m]=cvdcj[Mj(2*(i4d[0]-1)-1,j,kk,n)];
      else cvdcj[m]=.5*(cvdcj[Mj(ii-1,j,kk,n)]+cvdcj[Mj(ii+1,j,kk,n)]);
    }
    Loop(j,0,i4d[1]-1) Loop(k,0,i4d[2]-1) Loop(i,0,i4d[0])/* k interp i face */
    {
      ii=2*i; jj=2*j+1;
      m=Mk(ii,jj,k,n);
      if (cvdck[m]>=0) continue;
      if (i==0)
      { 
        if (csym[6]=='r' || csym[6]=='R') 
        { 
          cvdck[m]=.5*(cvdck[Mk(1,jj,k,n)]+cvdck[Mk(2*(i4d[0]-1)-1,jj,k,n)]);
          cvdck[Mk(2*(i4d[0]-1),jj,k,n)]=cvdck[m];
        }
        else cvdck[m]=cvdck[Mk(1,jj,k,n)];
      }
      else if (i==i4d[0]-1) cvdck[m]=cvdck[Mk(2*(i4d[0]-1)-1,jj,k,n)];
      else cvdck[m]=.5*(cvdck[Mk(ii-1,jj,k,n)]+cvdck[Mk(ii+1,jj,k,n)]);
    }
    
  }
}
/* ---------------------------- */
int cvdcfill(int *ipa, int *ipb, int iptc, int *i4d, int *i4dp, int iall, double *cvdc, 
             double *cvdci, double *cvdcj, double *cvdck, int *wherepx, int nyave)
{ 
  int ip[4],i,L=-1,ist[2]={0,0},iend[2]={0,0},iloop=0,ia,ib,warn=0, j;
  double cvdn[3];
  int bug=0;
  
  if (nyave>0) L=nyave-1; 
  else Loop(i,0,3) if (ipa[i]!=ipb[i]) L=i;
  
  Loop(i,0,4) ip[i]=ipa[i]; 
  
  if (bug>2) printout("normal","cvdcfill ipa %d %d %d to ipb %d %d %d L %d \n",
                    ipa[0],ipa[1],ipa[2],ipb[0],ipb[1],ipb[2],L); 
  cvdn[L]=.5;
  if (nyave==0)  /* one side of corner */
  {
    iend[0]=max(ipa[L],ipb[L]); 
    ist[0]=min(ipa[L],ipb[L]); 
    iloop=1;
    Loop(i,0,3) 
    if (i!=L)
    { 
      cvdn[i]=cvdc[iptc+i*iall]; 
      if (cvdn[i]<0) {cvdn[i]=.5; warn++;} 
    }
  }
  else /* average ends for linear fill */
  {  
    if (ipa[L]<ipb[L])
    { 
      ist[0]=ipa[L]; iend[0]=ipb[L]; iloop=1;
    }
    else /* crossing repeat */
    {
      ist[0]=0; iend[0]=ipb[L];
      ist[1]=ipa[L]; iend[1]=i4dp[L];
      iloop=2;
    }
    if (bug>2) printout("normal","iloop %d ist %d %d iend %d %d \n",
                      iloop,ist[0],ist[1],iend[0],iend[1]); 
    Loop(i,0,3) 
    if (i!=L)
    {  
      cvdn[i]=.5*(cvdc[In4low(ipa,i4d)+i*iall] +cvdc[In4low(ipb,i4d)+i*iall]);
      if (cvdn[i]<0) 
      {
        cvdn[i]=max(cvdc[In4low(ipa,i4d)+i*iall],cvdc[In4low(ipb,i4d)+i*iall]);
        warn++;
        if (cvdn[i]<0) cvdn[i]=.5;
      }
    }
  }
  if (bug>2) printout("normal","cvdn %lg %lg %lg warn %d\n",cvdn[0],cvdn[1],cvdn[2],warn);
  Loop(j,0,iloop)
  {
    Loop(ip[L],ist[j]+1,iend[j]) /* cvdc */
    {
      Loop(i,0,3) cvdc[In4low(ip,i4d)+i*iall]=cvdn[i];
      wherepx[In4(ip,i4dp)]=0;
      if (L==0) 
      {
        Loop(ia,0,3) Loop(ib,0,3) {cvdci[Mi(ip[0]-1,2*(ip[1]-1)+ia,2*(ip[2]-1)+ib,ip[3])]=cvdn[L]; 
          if (bug>2) printout("normal","Mi %d %d %d\n",ip[0]-1,2*(ip[1]-1)+ia,2*(ip[2]-1)+ib);  }
      }
      if (L==1) 
      {
        Loop(ia,0,3) Loop(ib,0,3) cvdcj[Mj(2*(ip[0]-1)+ia,ip[1]-1,2*(ip[2]-1)+ib,ip[3])]=cvdn[L]; 
      }
      if (L==2) 
      { 
        Loop(ia,0,3) Loop(ib,0,3) cvdck[Mk(2*(ip[0]-1)+ia,2*(ip[1]-1)+ib,ip[2]-1,ip[3])]=cvdn[L]; 
      }
    }
    if (L==0) 
    { 
      Loop(ip[0],2*ist[j],2*iend[j]-1)
      {
        cvdcj[Mj(ip[0],ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdn[1];
        cvdck[Mk(ip[0],2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdn[2];
        if (bug>2) printout("normal","Mj %d %d %d   Mk %d %d %d\n",
                          ip[0],ip[1]-1,2*(ip[2]-1)+1,ip[0],2*(ip[1]-1)+1,ip[2]-1); 
      }
    }
    else if (L==1)
    {
      Loop(ip[1],2*ist[j],2*iend[j]-1)
      { 
        cvdci[Mi(ip[0]-1,ip[1],2*(ip[2]-1)+1,ip[3])]=cvdn[0];
        cvdck[Mk(2*(ip[0]-1)+1,ip[1],ip[2]-1,ip[3])]=cvdn[2];
      }
    }
    else if (L==2)
    {
      Loop(ip[2],2*ist[j],2*iend[j]-1)
      {
        cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,ip[2],ip[3])]=cvdn[0];
        cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,ip[2],ip[3])]=cvdn[1];
      }
    }
    
  }
  return (warn);
}
/* --------------------------- */
/* fix corner interpolation */
void cvdccorner(int its, int ite, FILE *fprint)
{  
  int *i4d, *wherep; /* needed arays */
  double *cvdci,*cvdcj,*cvdck,*cvdc; 
  char *matchside; 
  
  int iall,iallp,ip[4],i4dp[4],ippt,is,js,ka,kd,i,ichange;
  double ci,cj,cit,cjt;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  wherep=(int *)need("wherep");
  cvdci=(double *)need("cvdcdouble");  
  cvdcj=cvdci+4*iall; 
  cvdck=cvdcj+4*iall;
  cvdc=(double *)need("cvdc");
  matchside=(char *)need("matchside");
  ichange=0;
  
  /* look for ij corners - non-parallel matching sides */
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)
  {
    ippt=In4(ip,i4dp);
    if (wherep[ippt]>=0) continue;
    if (matchside[ippt]=='j') {is=0; js=0;} /* -iside matches -jside */
    else if (matchside[ippt]=='J') {is=0; js=1;} /* -iside matches +jside */ 
    else if (matchside[ippt+iallp]=='j') {is=1; js=0;} /* +iside matches -jside */
    else if (matchside[ippt+iallp]=='J') {is=1; js=1;} /* +iside matches +jside */  
    else continue;
    for (ka=-2;ka<1;ka++)
    { kd=2*ip[2]+ka;
      ci=cvdci[Mi(ip[0]-1,(ip[1]-1+js)*2,kd,ip[3])];
      cj=cvdcj[Mj((ip[0]-1+is)*2,ip[1]-1,kd,ip[3])];
      if (is==js) cit=.5*(ci+cj);
      else cit=.5*(ci+(1-cj));
      if (cit!=ci) ichange++;
      if (is==js) cjt=cit; 
      else cjt=1-cit;
      cvdci[Mi(ip[0]-1,(ip[1]-1+js)*2,kd,ip[3])]=cit;
      cvdcj[Mj((ip[0]-1+is)*2,ip[1]-1,kd,ip[3])]=cjt;
      if (is==0) cit=1-sqrt(1-cit);
      else cit=sqrt(cit);
      if (js==0) cjt=1-sqrt(1-cjt);
      else cjt=sqrt(cjt);
      cvdci[Mi(ip[0]-1,(ip[1]-1)*2+1,kd,ip[3])]=cit;
      cvdcj[Mj((ip[0]-1)*2+1,ip[1]-1,kd,ip[3])]=cjt;
      cvdc[In4low(ip,i4d)]=cit;
      cvdc[In4low(ip,i4d)+iall]=cjt;
    }
  } 
    
}
/* --------------------------- */
/* set sleeping control volume centers and associated faces */
void cvdcsleep(int its, int ite, FILE *fprint, int itest)
{ 
  int *i4d,*wherep,*match,*matchpc; /* needed arays */
  double *cvdc,*cvdci,*cvdcj,*cvdck; char *cmatch;  
  int *wherepx; /* temporary array */
  
  int i,k=0,L=0,iall,iallp,ip[4],i4dp[4],ippt,ipa[2][4];
  int kk[4],iend[4],iptc,nycorner,icvd0,icvd1;
  int mat[2][2][2],ig[4],ia[3],m,nythrur;
  int two[3]={2,2,2};
  
  double cvda0,cvda1,cvda2;
  int  bug=0,warn=0;  
  int nyuser=1;
  char cx[9]="iIjJkKrR";
  char cxlin[7]="IiJjKk";
  
  i4d=(int *)need("idim4d");
  match=(int *)need("match");
  matchpc=(int *)need("matchpc");
  cmatch=(char *)need("matchside");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  cvdc=(double *)need("cvdc");
  cvdci=(double *)need("cvdcdouble");
  cvdcj=cvdci+4*iall; 
  cvdck=cvdcj+4*iall;
  wherep=(int *)need("wherep");
  wherepx=(int *)tmalloca(iallp,'i');
  Loop(i,0,iallp) wherepx[i]=wherep[i];
  
  /* do corners first with there sleeping extensions */
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)
  {
    ippt=In4(ip,i4dp);
    if (wherepx[ippt]>=0) continue;
    iptc=In4low(ip,i4d);
    nycorner=0;
    Loop(i,0,6) 
    {  
      if (cmatch[ippt+i*iallp]!='n' && cmatch[ippt+i*iallp]!=cxlin[i])
      { 
        kk[nycorner+1]=i; nycorner++; 
      }
    }
    if (nycorner==0) continue;
    nythrur=cmatchfollow(kk[1],cmatch,wherep,ippt,i4dp,iallp,(int)0,   &iend[0],&kk[0]);
    nythrur=cmatchfollow(kk[2],cmatch,wherep,ippt,i4dp,iallp,(int)0,   &iend[3],&kk[3]);
    if (iend[0]==0 || iend[3]==0) 
    {
      printout("warning cvdcsleep","apparent corner at %d %d %d is not a corner, ends %d %d\n",
              ip[0],ip[1],ip[2],iend[0],iend[3]);
      continue;
    }
    iexpand(iend[0],i4dp,ipa[0]);   
    iexpand(iend[3],i4dp,ipa[1]); 
    icvd0=In4low(ipa[0],i4d);   /* outside points */
    icvd1=In4low(ipa[1],i4d);
    L=kk[0]/2; k=kk[3]/2;
    i=(ipa[0][L]-ipa[1][L])*(ipa[0][k]-ipa[1][k]);
    cvda0=cvdc[icvd0+k*iall];
    cvda1=cvdc[icvd1+L*iall];
    if (i>0) 
    { 
      cvdc[iptc+L*iall]=.5*(1-cvda0+cvda1);
      cvdc[iptc+k*iall]=.5*(1-cvda1+cvda0);
    }
    else 
    {	
      cvdc[iptc+L*iall]=.5*(cvda0+cvda1);
      cvdc[iptc+k*iall]=.5*(cvda0+cvda1);
    }
    i=3-L-k;
    cvdc[iptc+i*iall]=.5*(cvdc[icvd0+i*iall]+cvdc[icvd1+i*iall]);
    /* copy i direction to 4 parallel faces */
    cvda2=cvdc[iptc+i*iall];
    Loop(m,0,3) ig[m]=ip[m]-1; ig[3]=ip[3];
    if (i==2)
    { 
      cvdck[Mk(2*ig[0],2*ig[1]+1,ig[2],ig[3])]=cvda2;
      cvdck[Mk(2*ig[0]+2,2*ig[1]+1,ig[2],ig[3])]=cvda2;
      cvdck[Mk(2*ig[0]+1,2*ig[1],ig[2],ig[3])]=cvda2;
      cvdck[Mk(2*ig[0]+1,2*ig[1]+2,ig[2],ig[3])]=cvda2;
    }
    
    if (bug>2)
    { 
      printout("normal","corner %d %d %d cvdc %lg %lg %lg ",
              ip[0],ip[1],ip[2],cvdc[iptc],cvdc[iptc+iall],cvdc[iptc+2*iall]);
      printout("normal","using %d %d %d cvdc %lg %lg %lg ",
				  ipa[0][0],ipa[0][1],ipa[0][2],cvdc[icvd0],cvdc[icvd0+iall],cvdc[icvd0+2*iall]);
      printout("normal","and %d %d %d cvdc %lg %lg %lg\n",
				  ipa[1][0],ipa[1][1],ipa[1][2],cvdc[icvd1],cvdc[icvd1+iall],cvdc[icvd1+2*iall]);
    }
    wherepx[ippt]=0;
    /* fill inside and faces */
    warn += cvdcfill(ipa[0],ip,iptc,i4d,i4dp,iall,cvdc,cvdci,cvdcj,cvdck,wherepx,0);
    warn += cvdcfill(ipa[1],ip,iptc,i4d,i4dp,iall,cvdc,cvdci,cvdcj,cvdck,wherepx,0);
    /* finish the corner, this part moved to cvdccorner  */	
   /* if (ip[L]>ipa[0][L]) cvdc[iptc+L*iall]=1.-sqrt(1.-cvdc[iptc+L*iall]);
    else cvdc[iptc+L*iall]=sqrt(cvdc[iptc+L*iall]);
    if (ip[k]>ipa[1][k]) cvdc[iptc+k*iall]=1.-sqrt(1.-cvdc[iptc+k*iall]);
    else cvdc[iptc+k*iall]=sqrt(cvdc[iptc+k*iall]);
    */
    cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc];
    cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc+iall];
    cvdck[Mk(2*(ip[0]-1)+1,2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdc[iptc+2*iall];
    
  }
  if (bug>0 || warn>0) printout("normal"," %d defaults of cvdc=.5 in resolving corners\n",warn);
  if (itest==-1) {printout("normal","corners only\n"); return;}
  
  /* fill in rest of sleeping volumes */	
  warn=0; 
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) 
  { 
    ippt=In4(ip,i4dp);
    if (wherepx[ippt]>=0) continue; /* already done */ 
    iptc=In4low(ip,i4d);
    nycorner=0;
    if (bug>2) 
      printout("normal","cvdc sleep pt at %d %d %d %d \n",ip[0]-1,ip[1]-1,ip[2]-1,ip[3]);
    
    kk[1]=-1;
    Loop(i,0,6) { if (cmatch[ippt+i*iallp]!='n') { kk[1]=i; break; }}
    if (kk[1]==-1) {printout("normal","didn't find side like should have at %d %d %d\n"
                            ,ip[0],ip[1],ip[2]); continue;}
    kk[2]=-1;
    Loop(k,0,6) { if (cmatch[ippt+k*iallp]==cx[i]) { kk[2]=k; break; }}
    if (kk[2]==-1) {printout("normal","didn't find side %c like should have at %d %d %d\n"
                            ,cx[i],ip[0],ip[1],ip[2]); continue;}
    L=kk[1]/2;
    nythrur=cmatchfollow(kk[1],cmatch,wherep,ippt,i4dp,iallp,nyuser,   &iend[0],&kk[0]);
    nythrur=cmatchfollow(kk[2],cmatch,wherep,ippt,i4dp,iallp,nyuser,   &iend[3],&kk[3]);
    if (bug>2) 
    { 
      printout("normal"," sleep at %d %d %d, matchside",ip[0],ip[1],ip[2]);
      printout("normal"," kk1 %d kk2 %d iend0 %d iend3 %d\n",kk[1],kk[2],iend[0],iend[3]);
    }
    if (iend[0]>0) iexpand(iend[0],i4dp,ipa[0]);
    if (iend[3]>0) iexpand(iend[3],i4dp,ipa[1]); 
    if (iend[0]>0 && iend[3]>0)  /* two points */
      warn += cvdcfill(ipa[0],ipa[1],iptc,i4d,i4dp,iall,cvdc,cvdci,cvdcj,cvdck,wherepx,L+1);
    else if (iend[0]<=0 && iend[3]<=0) continue;
    else if (iend[0]>0) 
      warn += cvdcfill(ipa[0],ip,In4low(ipa[0],i4d),i4d,i4dp,iall,cvdc,cvdci,cvdcj,cvdck,wherepx,0);
    else 
      warn += cvdcfill(ipa[1],ip,In4low(ipa[1],i4d),i4d,i4dp,iall,cvdc,cvdci,cvdcj,cvdck,wherepx,0);
  } /* ip loop */
  if (bug>0 || warn>0) printout("warning cvdcsleep"," %d defaults of cvdc=.5 in resolving line sleeps\n",warn);
  
  /* if not finished use matchpc to finish */
  warn=0; k=0;
  Loop(m,0,10)
  { 
    warn=0; k=0;
    Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)
    {
      ippt=In4(ip,i4dp);
		if (wherepx[ippt]>=0) continue;
		warn++;
		iexpand(matchpc[ippt],i4dp,ipa[0]);
		iexpand(matchpc[ippt+iallp],i4dp,ipa[1]);
		icvd0=In4low(ipa[0],i4d);
		icvd1=In4low(ipa[1],i4d);
		if (cvdc[icvd0]<0 || cvdc[icvd1]<0) continue;
		warn--; k++; wherepx[ippt]=0;
		if (bug>2) printout("normal","matchpc-resolve of %d %d %d using %d %d %d and %d %d %d\n",
                         ip[0],ip[1],ip[2],ipa[0][0],ipa[0][1],ipa[0][2],ipa[1][0],ipa[1][1],ipa[1][2]); 
		iptc=In4low(ip,i4d);
		if (wherep[matchpc[ippt]]<0) Loop(L,0,3) cvdc[iptc+L*iall]=cvdc[icvd0+L*iall];
		else if (wherep[matchpc[ippt+iallp]]<0) Loop(L,0,3) cvdc[iptc+L*iall]=cvdc[icvd1+L*iall];
		else Loop(L,0,3) cvdc[iptc+L*iall]=.5*(cvdc[icvd0+L*iall]+cvdc[icvd0+L*iall]);
		cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc];
		cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc+iall];
		cvdck[Mk(2*(ip[0]-1)+1,2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdc[iptc+2*iall];
    }
    if (bug>0) printout("normal"," %d resolved with matchpc, %d still unresolved\n",k,warn);
    if (warn==0 || k==0) break;
  }
  warn=0; /* count to see if finish ok  */
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)   /* fill in rest of sleeping control volumes */
  {
    ippt=In4(ip,i4dp);
    if (wherepx[ippt]>=0) continue; /* already done */
    printout("normal"," not resolved at %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]);
    warn++;
    iptc=In4low(ip,i4d);
    Loop(L,0,3) cvdc[iptc+L*iall]=.5;
    cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc];
    cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdc[iptc+iall];
    cvdck[Mk(2*(ip[0]-1)+1,2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdc[iptc+2*iall];
  }
  if (bug>0 || warn>0) printout("normal"," %d defaults of cvdc=.5 for sleeping c.v.\n",warn);
  if (itest==-2) {printout("normal","corners and cmatch sleeps only\n"); return;}
  warn=0;
  /* check for missed sleeping i, j and k faces */
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)
  { 
    if (wherep[In4(ip,i4dp)]>=0) continue;
    ig[3]=ip[3];
    Loop(i,0,3) ig[i]=ip[i]-1;
    cvda0=cvdci[Mi(ig[0],2*ig[1]+1,2*ig[2]+1,ig[3])];
    cvda1=cvdcj[Mj(2*ig[0]+1,ig[1],2*ig[2]+1,ig[3])];
    cvda2=cvdck[Mk(2*ig[0]+1,2*ig[1]+1,ig[2],ig[3])];
    if (bug>2) printout("normal","check face %d %d %d, cvdcenter %lg %lg %lg mat",
                      ip[0],ip[1],ip[2],cvda0,cvda1,cvda2); 
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
      mat[ia[0]][ia[1]][ia[2]]=match[In4(ig,i4d)];
      if (bug>2) printout("normal"," %d",mat[ia[0]][ia[1]][ia[2]]); 
    }
    Loop(i,0,3) ig[i]=ip[i]-1;
    if (mat[0][0][0]==mat[1][0][0] && mat[0][1][0]==mat[1][1][0] 
        && mat[0][0][1]==mat[1][0][1] && mat[0][1][1]==mat[1][1][1])
    { 
      /* iface match  copy j and k from center */
      m=Mj(2*ig[0],ig[1],2*ig[2]+1,ig[3]); if (cvdcj[m]<0) {cvdcj[m]=cvda1; warn++;}
      m=Mj(2*ig[0]+2,ig[1],2*ig[2]+1,ig[3]); if (cvdcj[m]<0) {cvdcj[m]=cvda1; warn++;}
      m=Mk(2*ig[0],2*ig[1]+1,ig[2],ig[3]); if (cvdck[m]<0) {cvdck[m]=cvda2; warn++;}
      m=Mk(2*ig[0]+2,2*ig[1]+1,ig[2],ig[3]); if (cvdck[m]<0) {cvdck[m]=cvda2; warn++;}
    }
    else if (bug>2) printout("normal"," not i"); 
    if (mat[0][0][0]==mat[0][1][0] && mat[1][0][0]==mat[1][1][0] 
        && mat[0][0][1]==mat[0][1][1] && mat[1][0][1]==mat[1][1][1])
    { 
      /* jface match  copy i and k from center */
      m=Mi(ig[0],2*ig[1],2*ig[2]+1,ig[3]); if (cvdci[m]<0) {cvdci[m]=cvda0; warn++;}
      m=Mi(ig[0],2*ig[1]+2,2*ig[2]+1,ig[3]); if (cvdci[m]<0) {cvdci[m]=cvda0; warn++;}
      m=Mk(2*ig[0]+1,2*ig[1],ig[2],ig[3]); if (cvdck[m]<0) {cvdck[m]=cvda2; warn++;}
      m=Mk(2*ig[0]+1,2*ig[1]+2,ig[2],ig[3]); if (cvdck[m]<0) {cvdck[m]=cvda2; warn++;}
    }
    else if (bug>2) printout("normal"," not j"); 
    if (mat[0][0][0]==mat[0][0][1] && mat[1][0][0]==mat[1][0][1] 
        && mat[0][1][0]==mat[0][1][1] && mat[1][1][0]==mat[1][1][1])
    { 
      /* kface match copy i and j from center */
      m=Mi(ig[0],2*ig[1]+1,2*ig[2],ig[3]); if (cvdci[m]<0) {cvdci[m]=cvda0; warn++;}
      m=Mi(ig[0],2*ig[1]+1,2*ig[2]+2,ig[3]); if (cvdci[m]<0) {cvdci[m]=cvda0; warn++;}
      m=Mj(2*ig[0]+1,ig[1],2*ig[2],ig[3]); if (cvdcj[m]<0) {cvdcj[m]=cvda1; warn++;}
      m=Mj(2*ig[0]+1,ig[1],2*ig[2]+2,ig[3]); if (cvdcj[m]<0) {cvdcj[m]=cvda1; warn++;}
    }
    else if (bug>2) printout("normal"," not k");
    if (bug>2) printout("normal","\n"); 
  } 
  if (bug>0 || warn>0) printout("normal"," set %d missed sleeping faces\n",warn); 
}
/* ------------------------------- */
/* read parameters for control volume limits */
void c_cvdcparm(FILE *fpin, FILE *fprint)
{
  int i;
  if (flimd>0) {free(flimd); free(climd);}
  flim=readdouble(fpin);   if (flim>0.5) flim=0.5;
  limd=readint(fpin);
  printout("normal","cvdc>= %lg and cvdc<= %lg \n",flim,1.-flim);
  if (limd>0)
  {
    flimd=(double*)smalloca(limd,'d');
    climd=(char*)smalloca(limd*2,'c');
    Loop(i,0,limd)
    {
      climd[i]=read1charname(fpin);
      climd[i+limd]=read1charname(fpin);
      flimd[i]=readdouble(fpin);
      printout("normal","  from %c to %c cvdc>= %lg\n",climd[i],climd[i+limd],flimd[i]);
    }
  }
}
/* ------------------------------------- */
/* read fixed parameters and initialize cvd centered with pt type limiitations */
void c_cvdcinit(FILE *fpin, FILE *fprint)
{
  int *i4d,*wherep; char *clt; /* needed arrays */
  double *cvdc, *cvdcdouble;  /* create and initialize */
  double *cvdci,*cvdcj,*cvdck;
  int i,iall,ip[4],ipt,ipr[4],i4dp[4];
  double cvdn[3];
  int itest=0;  /* set to zero when debugged */
  
  c_cvdcparm(fpin,fprint);
  /*printout("normal"," itest %d\n",itest); */
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  clt=(char *)need("clt");
  ipr[0]=1; Loop(i,1,4) ipr[i]=ipr[i-1]*i4d[i-1];
  iall=ipr[3];
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  
  cvdc=(double *)createarray("cvdc",iall*3,'d',0); /* create and initialize as -1 so can see if set */
  Loop(i,0,iall*3) cvdc[i]=-1;
  cvdcdouble=(double *)createarray("cvdcdouble",iall*12,'d',0);
  Loop(i,0,iall*12) cvdcdouble[i]=-1;
  cvdci=cvdcdouble; cvdcj=cvdci+4*iall; cvdck=cvdcj+4*iall;
  
  /* set to .5 with limits */
  Loop(ip[3],0,i4d[3]) Loop3(ip,1,i4d) 
  { 
    Loop(i,0,3)  cvdn[i]=.5;
    ipt=In4low(ip,i4d);
    cvdclimit1(cvdn,ipt,ipr,clt);
    Loop(i,0,3) cvdc[ipt+iall*i]=cvdn[i];
    cvdci[Mi(ip[0]-1,2*(ip[1]-1)+1,2*(ip[2]-1)+1,ip[3])]=cvdn[0];
    cvdcj[Mj(2*(ip[0]-1)+1,ip[1]-1,2*(ip[2]-1)+1,ip[3])]=cvdn[1];
    cvdck[Mk(2*(ip[0]-1)+1,2*(ip[1]-1)+1,ip[2]-1,ip[3])]=cvdn[2];
  }
  if (itest==1) return;
  cvdcsleep(0,i4d[3]-1,fprint,itest);
  if (itest==2 || itest==-1) return;
  cvdcface(0,i4d[3]-1,fprint);
  if (itest==3) return;
  cvdclines(0,i4d[3]-1,fprint);
  cvdccorner(0,i4d[3]-1,fprint);
}
/* ------------------------- */
/* reset control volumes based on convection, and possibly viscosity and time */
void c_cvdcreset(FILE *fpin, FILE *fprint)
{
  int *i4d; double *cvdcdouble;
  char *nameisovisc, *namedt; double t2cvlim; /* input parameters */
  int *itrange; double *visc, *dt;/* use if available */
  int its,ite,i,iall,sizedt=0;
  int itest=0;
  
  nameisovisc=readname(fpin);  /* input */
  namedt=readname(fpin);
  t2cvlim=readdouble(fpin);
  /*itest=readint(fpin);  */ /* comment out when debugged */
  /*printout("normal","itest %d\n",itest); */
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  cvdcdouble=(double *)need("cvdcdouble");
  Loop(i,0,iall*12) cvdcdouble[i]=-1;
  
  printout("normal","viscosity %s, timestep %s, t2cvlim %lg",nameisovisc,namedt,t2cvlim);
  visc=(double *)find(nameisovisc);
  dt=(double *)find(namedt);
  if (visc>0) {printout("normal","  using viscosity"); }
  if (dt>0) 
  {
    sizedt=arraysize(namedt); 
    printout("normal","  using timestep %d",sizedt); 
  }
  printout("normal","\n");
  
  cvdcactive(its,ite,visc,dt,sizedt,t2cvlim,fprint);
  /*printout("normal","active done\n"); */
  if (itest==1) return;
  cvdcsleep(its,ite,fprint,itest);
  printout("normal","sleep done\n");
  if (itest==2 || itest==-1 || itest==-2) return;
  cvdcface(its,ite,fprint);
  printout("normal","face done\n");
  if (itest==3) return;
  cvdclines(its,ite,fprint);
  cvdccorner(its,ite,fprint);  /* check corners */
}
