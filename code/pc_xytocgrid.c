/* contains pc_xytocgrid */
#include "global.h"
#include "p_plot.h"
/*  c-grid around the blade between the points 
 input: append, 
 cs, ce, numasbe, (as, ae, bs, be, ) numasbe times
 distance, nptsn, expand (neg for flat + for c grid),
 prop until prop="",
 */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Int2d(p,n,id) ( (1-fin[0])*(1-fin[1])*(p)[n+ilowin[0]+ilowin[1]*id[0]]+fin[0]*(1-fin[1])*(p)[n+ilowin[0]+1+ilowin[1]*id[0]]+(1-fin[0])*fin[1]*(p)[n+ilowin[0]+(ilowin[1]+1)*id[0]]+fin[0]*fin[1]*(p)[n+ilowin[0]+1+(ilowin[1]+1)*id[0]] )
/* ---------------------------- */
void pc_xytocgrid(FILE *fpin, FILE *fprint)
{
  int *i4d; double *a[4],*x[3],*u[3]; char *clt, *csym; /* need */
  char *append; double cs,ce; int numab; double as[20],ae[20],bs[20],be[20];
  double dist; int nptsn; double expand;
  char *name;   /* input parameters */
  double *v;
  int *i4dnew; double *anew[4], *xnew[3],*unew[3]; char *cltnew=0; double *var; /* create */
  double *xm[2];
  int iall,i,j,k,is[20],ie[20],js[20],je[20],npts,istr,inormmax=0,signs,signn;
  int ks,ke,kmid,ii[4],nyok,n,m,mm,iorj,iallnew,iin[4];
  char namg[7][20]={"idim4d","abcd","xyz","clt","U1","U2","U3"},namenew[50];
  double f,ds,*dxy,*dxyn,*fin,*fout,dn,us,un,ff;
  int *ilowin, *ilowout;
   Array *array;
  int *inew; double *fnew, *dists, *distn, *vecstr, *vecnorm; /* temporary arrays */
  int bug=0;
  
  /* get needed arrays */
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double *)need("abcd");
  Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[1]+iall;
  clt=(char *)need("clt");
  csym=(char *)need("csym");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  
  /* read input except for prop */
  
  append=readname(fpin);
  cs=readdouble(fpin);
  ce=readdouble(fpin);
  ks=findex(cs,a[2],i4d[2],&f)+f+.5;
  ke=findex(ce,a[2],i4d[2],&f)+f+.5;
  kmid=(ks+ke)/2;
  numab=readint(fpin);
  numab=min(numab,20);
  printout("normal"," append %s, c=%lg thru %lg, (k= %d ,%d) numab=%d\n",append,cs,ce,ks,ke,numab);
  npts=0;
  for (i=0;i<numab;i++) 
  {
    as[i]=readdouble(fpin);
    ae[i]=readdouble(fpin);
    bs[i]=readdouble(fpin);
    be[i]=readdouble(fpin);
    is[i]=findex(as[i],a[0],i4d[0],&f)+f+.5;
    ie[i]=findex(ae[i],a[0],i4d[0],&f)+f+.5;
    js[i]=findex(bs[i],a[1],i4d[1],&f)+f+.5;
    je[i]=findex(be[i],a[1],i4d[1],&f)+f+.5;
    
    printout("normal","a = %lg, %lg,  b = %lg, %lg gives i=%d, %d, j=%d,%d\n",
            as[i],ae[i],bs[i],be[i],is[i],ie[i],js[i],je[i]);
    npts+=abs(ie[i]-is[i])+abs(je[i]-js[i]);
  }
  dist=readdouble(fpin);
  nptsn=readint(fpin);
  expand=readdouble(fpin);
  printout("normal","norm dist= %lg, maxnormal=%d, expandsion factor %lg\n",dist,nptsn,expand);
  printout("normal","npts around surface= %d\n",npts);
  
  /* temporary arrays */
  inew=(int *)tmalloca(2*npts*nptsn,'i');
  fnew=(double *)tmalloca(2*npts*nptsn,'d');
  dists=(double *)tmalloca(npts,'d');
  vecstr=(double *)tmalloca(2*npts,'d');
  vecnorm=(double *)tmalloca(2*npts,'d');
  distn=(double *)tmalloca(npts*nptsn,'d');
  
  /* loop streamwise direction sides */
  istr=0;
  ii[2]=kmid; ii[3]=0;
  xm[0]=x[0]+kmid*i4d[0]*i4d[1];
  xm[1]=x[1]+kmid*i4d[0]*i4d[1];
  
  ds=0;
  for (n=0;n<numab;n++)
  {
    ii[0]=is[n]; ii[1]=js[n];
    if (ie[n]>is[n]) {signs=1; signn=1; iorj=0; k=ie[n]-is[n]; }
    else if (ie[n]<is[n]) {signs=-1; signn=-1; iorj=0; k=is[n]-ie[n]; }
    else if (je[n]>js[n]) {signs=1; signn=-1; iorj=1; k=je[n]-js[n]; }
    else if (je[n]<js[n]) {signs=-1; signn=1; iorj=1; k=js[n]-je[n]; }
    else
    {
      printout("warning pc_xytocgrid","single point ignored\n");
      continue;
    }
    /* loop over points along side */
    for (i=0;i<k;i++) 
    {
      dists[istr]=0;
      if (istr>0) dists[istr]=dists[istr-1];
      /* addresses */
      dxy=vecstr+2*istr;
      dxyn=vecnorm+2*istr;
      ilowin=inew+2*istr;
      fin=fnew+2*istr;
      
      for (j=0;j<2;j++)
      {
        ii[iorj]+=signs;
        dxy[j]=x[j][In4(ii,i4d)];
        ii[iorj]-=signs;
        dxy[j]-=x[j][In4(ii,i4d)];
      }
      dists[istr]+=.5*ds;
      ds=sqrt(dxy[0]*dxy[0]+dxy[1]*dxy[1]);
      if (bug>0) printout("normal","istr %d dxy str %lg %lg, ds %lg,",istr,dxy[0],dxy[1],ds);
      dists[istr]+=.5*ds;
      dxy[0]/=ds; dxy[1]/=ds;
      dxyn[0]=-dxy[1];
      dxyn[1]=dxy[0];
      if (bug>0) printout("normal"," dxy norm %lg %lg\n",dxyn[0],dxyn[1]);
      ilowin[iorj]=ii[iorj]; if (signs<0) ilowin[iorj]--;
      fin[iorj]=.5;
      ilowin[1-iorj]=ii[1-iorj]; fin[1-iorj]=0;
      if (signn<0) {ilowin[1-iorj]--; fin[1-iorj]=1; }
      distn[istr]=0;
      /* loop normal to wall */
      for (m=1;m<nptsn;m++)  
      {
        /* addresses */
        ilowout=inew+2*istr+2*npts*m;
        fout=fnew+2*istr+2*npts*m;
        
        for (j=0;j<2;j++) {ilowout[j]=ilowin[j]; fout[j]=fin[j]; }
        nyok=1;
        dn=p_nextpt2d(dxyn,xm,i4d,ilowin,fin,ilowout,fout); 
        /* do checks new dn, same, bounds */
        if (dn<0) nyok=0;
        if (fin[0]==fout[0] && ilowin[0]==ilowout[0]
            && fin[1]==fout[1] && ilowin[1]==ilowout[1]) nyok=0;
        if (ilowout[0]<0 || ilowout[0]> i4d[0]-2) nyok=0;
        if (ilowout[1]<0 || ilowout[1]> i4d[1]-2) nyok=0;
        distn[istr+m*npts]=distn[istr+(m-1)*npts]+dn;
        if ( nyok==1 && distn[istr+m*npts]>=dist)
        {
          f=(dist-distn[istr+(m-1)*npts])/dn;
          for (j=0;j<2;j++) 
          {
            ff=(1-f)*(fin[j]+ilowin[j])+f*(fout[j]+ilowout[j]);
            ilowout[j]=ff;
            if (ilowout[j]==i4d[j]-1) ilowout[j]--;
            fout[j]=ff-ilowout[j];
          }
          distn[istr+m*npts]=dist;
          m++;
          fin=fout;
          ilowin=ilowout;
          nyok=-1;
        }
        if (nyok==1 && dn==0)
        {
          fin[0]=fout[0]; fin[1]=fout[1];
          ilowin[0]=ilowout[0]; ilowin[1]=ilowout[1];
          if (bug>0) printout("normal"," ok dn=0 %d %lg %d %lg\n",ilowin[0],fin[0],ilowin[1],fin[1]);
          m--;
        }
        else if (nyok==1)
        {
          inormmax=max(inormmax,m);
          fin=fout;
          ilowin=ilowout;
          if (bug>0) printout("normal"," ok %d %lg %d %lg\n",ilowin[0],fin[0],ilowin[1],fin[1]);
        }
        else
        {
          
          if (bug>0) printout("normal"," end m loop m %d max %d nyok %d ilow f %d %lg %d %lg\n",m,inormmax,nyok,
                 ilowout[0],fout[0],ilowout[1],fout[1]);
          inormmax=max(inormmax,m);
          for (mm=m;mm<nptsn;mm++)
          { 
            distn[istr+mm*npts]=distn[istr+(m-1)*npts];
            ilowout=inew+2*istr+2*npts*mm;
            fout=fnew+2*istr+2*npts*mm;
            for (j=0;j<2;j++) {ilowout[j]=ilowin[j]; fout[j]=fin[j]; }
          }
          m=nptsn+1;
          
        } 
      } /* end normal to wall */
      istr++;
      ii[iorj]+=signs;
    }
  }
  printout("normal","inormmax %d\n",inormmax);
  /* set up new grid  and velocity components */
  Loop(i,0,7) strncat(namg[i],append,strlen(append)+1);
  i4dnew=(int *)createarray(namg[0],4,'i',1);
  i4dnew[0]=npts;
  i4dnew[1]=inormmax;
  i4dnew[2]=ke-ks+1;
  i4dnew[3]=1;
  iallnew=Prod4(i4dnew);
  /* set abcdnew */
  anew[0]=(double *)createarray(namg[1],i4dnew[0]+i4dnew[1]+i4dnew[2]+i4dnew[3],'d',1);
  Loop(i,1,4) anew[i]=anew[i-1]+i4dnew[i-1];
  Loop(i,0,npts) anew[0][i]=dists[i];
  Loop(i,0,i4dnew[1]) anew[1][i]=i*dist/i4dnew[1];
  Loop(i,0,i4dnew[2]) anew[2][i]=a[2][i+ks];
  anew[3][0]=0;
  xnew[0]=(double *)createarray(namg[2],iallnew*3,'d',1);
  xnew[1]=xnew[0]+iallnew; xnew[2]=xnew[1]+iallnew;
  cltnew=(char *)createarray(namg[3],iallnew,'c',1);
  Loop(i,0,3) unew[i]=(double *)createarray(namg[4+i],iallnew,'d',1);
  /* set new grid and velocities */
  iin[3]=0;
  for (iin[0]=0;iin[0]<i4dnew[0];iin[0]++)
    for (iin[1]=0;iin[1]<i4dnew[1];iin[1]++)
    { 
      ilowin=inew+2*iin[0]+2*npts*iin[1];
      fin=fnew+2*iin[0]+2*npts*iin[1];
      for (iin[2]=0;iin[2]<i4dnew[2];iin[2]++)
      {
        mm=In4(iin,i4dnew);
        if (iin[1]==0) cltnew[mm]='w'; else cltnew[mm]='f';
        n=i4d[0]*i4d[1]*(iin[2]+ks);
        if (bug>0) printout("normal","iin %d %d %d ilow %d %d\n",iin[0],iin[1],iin[2],ilowin[0],ilowin[1]);
        for (j=0;j<3;j++) xnew[j][mm]=Int2d(x[j],n,i4d);
        for (j=0;j<3;j++) unew[j][mm]=Int2d(u[j],n,i4d);
        if (expand<0) /* s-n grid */
        {
          xnew[0][mm]=dists[iin[0]];
          xnew[1][mm]=distn[iin[0]+npts*iin[1]];
          us=unew[0][mm]*vecstr[2*iin[0]]+unew[1][mm]*vecstr[2*iin[0]+1];
          un=unew[0][mm]*vecnorm[2*iin[0]]+unew[1][mm]*vecnorm[2*iin[0]+1];
          unew[0][mm]=us;
          unew[1][mm]=un;
        }
      }
    }
  /* extra properties */
  printout("normal"," converting arrays:\n");
  n=100;
  while (n)
  { 
    name=readname(fpin);
    if (name[0]=='\0') break;
    if (strcmp(name,"c:")==0) break;
    n--;
    array=findarray(name);
    if (array==0)  {printout("warning"," array %s not found,\n",name); continue; }
   
    if (array->type!='d' || (array->size)%iall !=0) 
    { printout("warning"," not on points type d, %s ignored\n",name); continue; }
    printout("normal"," %s",name);
    j=(array->size)/iall;
    v=(double *)array->pointer;
    strcpy(namenew,name);
    strncat(namenew,append,strlen(append)+1);
    var=(double *)createarray(namenew,j*iallnew,'d',1);
    printout("normal","  %s\n",namenew);
    for (iin[0]=0;iin[0]<i4dnew[0];iin[0]++)
      for (iin[1]=0;iin[1]<i4dnew[1];iin[1]++)
      { 
        ilowin=inew+2*iin[0]+2*npts*iin[1];
        fin=fnew+2*iin[0]+2*npts*iin[1];
        for (iin[2]=0;iin[2]<i4dnew[2];iin[2]++)
        {
          mm=In4(iin,i4dnew);
          if (iin[1]==0) cltnew[mm]='w'; else cltnew[mm]='f';
          n=i4d[0]*i4d[1]*(iin[2]+ks);
          var[mm]=Int2d(v,n,i4d);
        } 
      }
  }
}