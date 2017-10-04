/* contains c_momrhsr */
#include "global.h"
/* rhs of momentum equations
 nynew        nynew 0,1 to start with current rhs or clean rhs
 read jcoef   integer, coefficient set to be used for rhs, if less than zero it is omitted 
 then one or more of the following names (in any order) for terms to be included, ending with ""
 bij          to include Reynolds stress terms, will need bij and qturb to be included
 dpdx         for overall dpdx, will need dpdx to be included
 zrotation    for coordinate rotation needs zrotation to be included and velocity
 vname        name starting with v for the viscous term d/dxi(vname dUi/dxj), include only one
 pname        name starting with p for the pressure dP/dxj term
 may include one p-pts pressure and one g-pts pressure
 terms are evaluated as surface integrals except rotation which is a volume integral  
 */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_momrhsr(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*noindfwpts;  /* needed arrays */
  double *rhs[3]={0,0,0};  /* set rhsU1 rhsU2 rhsU3 */
  int *itrange;  /* use if available */
  
  int *cpflop_n, **cpflop_i, *wherefw=0, *wherep=0, *nocoefs; /* arrays that might be needed */
  double **cpflop_c, *rho=0, *u[3]={0,0,0}, *dudxi[3]={0,0,0}, *x[3]={0,0,0}, *q=0,*qmid=0;	
  char *csym=0; 
  double *dpdx=0, *bij=0, *pp=0, *pg=0, *visc=0, *zrotation,*rhomid=0, *bijm=0; /* use if wanted or available */
  
  double fac=1; /* fixed parameter */
  int jcoef=0; /* read in this version */
  char *name,oldnew;
  
  double zrot=0,xrot=0, yrot=0; int nycentrifugal=0;/* zrotation, about xrot,yrot */
  
  int i,j,m,iw[8],its,ite,iall,i4dp[4],i4dm[4],iallp,iallm,mmid;
  int ip[4],id[4],ig[4],ia[3],L,L1,L2,iim[4],iip[4],iptm,iptp,ipt[8],iwm,iwp;
  int n=100;
  int n2ii[6] = {0,1,2,0,0,1};
  int n2jj[6] = {0,1,2,1,2,2};
  int two[3] = {2,2,2};
  double area[3],xmid[3],qm,rhom,uiuj[3][3],drhs;
  double vol[8],fmid[8][3],xm[2],rum[2],ff[3],ffc[8];
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  Loop(i,0,3) {i4dp[i]=i4d[i]+1; i4dm[i]=i4d[i]-1; }
  i4dp[3]=i4d[3]; i4dm[3]=i4d[3];
  iallp=Prod4(i4dp);  iallm=Prod4(i4dm);
  noindfwpts=(int *)need("noindfwpts");
  
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  
  oldnew=read1charname(fpin);   /* read input */
  
  jcoef=readint(fpin); 
  printout("normal","oldnew %c, jcoef %d",oldnew,jcoef);
  if (jcoef>=0) 
  { 
    nocoefs=(int *)need("nocoefs");
    if (jcoef<nocoefs[0]) printout("normal"," included,");
    else {printout("normal"," omitted,"); jcoef=-1;}
  }
  while (n)   /* read names of what to include */
  {
    name=readname(fpin);
    if (name[0]=='\0')  break;
    if (strcmp(name,"c:")==0) { fseek(fpin,(long)(-3),1); break; }
    n--;	
    printout("normal"," %s ",name);
    if (strcmp(name,"dpdx")==0)
    {  
      dpdx=(double *)find("dpdx");
      if (dpdx>0 && dpdx[0]!=0) printout("normal","=%lg included,",dpdx[0]);
      else {dpdx=0; printout("normal"," omitted,"); }
    }
    else if (strcmp(name,"bij")==0)
    {	
      if (arraysize("qturb")==iall) q=(double *)need("qturb");
      if (arraysize("bij")==6*iall) bij=(double *)need("bij");
      if (bij>0 && q>0) printout("normal"," included,");
      else {bij=0; printout("normal"," omitted,");}
    }
    else if (strcmp(name,"bijm")==0)
    {	
      if (arraysize("qmid")==iallm) qmid=(double *)need("qmid");
      if (arraysize("bijm")==6*iallm) bijm=(double *)need("bijm");
      if (arraysize("rhom")==iallm) rhomid=(double *)need("rhom");
      if (bijm>0 && qmid>0 && rhomid>0) printout("normal"," included,");
      else {bijm=0; printout("normal"," omitted,");}
    }
    else if (strcmp(name,"zrotation")==0)
    { 
      zrotation=(double *)find("zrotation");
      if (zrotation>0 && zrotation[0]!=0)
      { 
        zrot=zrotation[0]; 
        printout("normal","=%lg included,",zrot); 
        if (arraysize("zrotation")>=3) 
        { 
          xrot=zrotation[1]; yrot=zrotation[2]; 
          printout("normal","about x=%lg, y=%lg,",xrot,yrot);
          nycentrifugal=1;
        }
      }
      else printout("normal"," omitted,");
    }
    else if (name[0]=='p')
    { 
      j=arraysize(name);
      if (j==iallp) { pp=(double *)need(name); printout("normal"," at ppts included,");}
      else if (j==iall) {pg=(double *)need(name); printout("warning c_momcamr"," at gpts not yet coded,");}
      else printout("normal"," omitted,");
    }
    else if (name[0]=='v')
    {
      if (arraysize(name)==iallm) {visc=(double *)need(name);}
		else printout("normal"," omitted,");
    }
  }
  printout("normal","\n");
  if (jcoef>=0) /* convection etc */
  {
    coefrhs(fprint,"U1","rhsU1",jcoef,oldnew,fac);
    coefrhs(fprint,"U2","rhsU2",jcoef,oldnew,fac);
    coefrhs(fprint,"U3","rhsU3",jcoef,oldnew,fac);
  }
  if (jcoef>=0 || oldnew=='o')
  { 
    rhs[0]=(double *)find("rhsU1");
    rhs[1]=(double *)find("rhsU2");
    rhs[2]=(double *)find("rhsU3");
  }
  if (rhs[0]==0) 
  {	
    rhs[0]=(double *)createarray("rhsU1",noindfwpts[0],'d',0);
    rhs[1]=(double *)createarray("rhsU2",noindfwpts[0],'d',0);
    rhs[2]=(double *)createarray("rhsU3",noindfwpts[0],'d',0);
    Loop(i,0,noindfwpts[0]) Loop(j,0,3) rhs[j][i]=0;
  }
  /* pressure at p-pts */
  if (pp>0)
  {
    cpflop_n=(int*)need("cpflop_n"); 
    cpflop_i=(int**)need("cpflop_i"); 
    cpflop_c=(double**)need("cpflop_c"); 
    Loop(iwp,noindfwpts[its+1],noindfwpts[ite+2])   
    {
      n=cpflop_n[iwp];
      if (n==0) continue;
      Loop(i,0,n) Loop(j,0,3) rhs[j][iwp]-=cpflop_c[iwp][i+n*j]*pp[cpflop_i[iwp][i]];
    }
  }
  /* gather other needed arrays */
  if (dpdx>0 || bij>0 || zrot!=0 || pg>0 || visc>0 || bijm>0)
  { 
    wherep=(int *)need("wherep");
    wherefw=(int *)need("wherefw");
    rho=(double *)need("rho");
    csym=(char *)need("csym");
    geom8init();
    
    if (visc>0)
    { 
      dudxi[0]=(double *)find("U1ddxi");
      dudxi[1]=(double *)find("U2ddxi");
      dudxi[2]=(double *)find("U3ddxi");
      if (dudxi[0]==0 || dudxi[1]==0 || dudxi[2]==0)
      {
        printout("warning c_momrhsr"," viscous term omitted, U1ddxi, U2ddxi, or U3ddxi, do not exist\n");
        visc=0;
      }
    }
    if (zrot!=0)
    { 
      u[0]=(double *)need("U1");
      u[1]=(double *)need("U2");
      u[2]=(double *)need("U3");
      x[0]=(double *)need("xyz");
      x[1]=x[0]+iall; x[2]=x[1]+iall;
    } 
  }
  
  /* surface integral sources */
  if (dpdx>0 || bij>0  || pg>0 || visc>0 || bijm>0)
  {
    Loop(ip[3],its,ite+1) 
    Loop(L,0,3)  /* surface(s) in each direction */
    {  
      L1=(L+1)%3; L2=(L+2)%3; 
      Loop(ip[L1],1,i4d[L1]) Loop(ip[L2],1,i4d[L2]) /* each c.v. */
      Loop(ip[L],0,i4dp[L])
      {
        id[3]=ip[3]; iim[3]=ip[3]; iip[3]=ip[3];
        iim[L]=ip[L]-1; iip[L]=ip[L];
        ia[L]=1;
        if (ip[L]==0) ia[L]=2;   /* boundaries */
        else if (ip[L]==i4d[L]) ia[L]=0;
        else if (wherep[In4(ip,i4dp)]<0) continue; /* not valid cont.c.v. */
        else if (csym[6+2*L]=='R') continue; /* omit for 2d repeat, except boundaries */
        
        Loop(ia[L1],0,2) Loop(ia[L2],0,2) /*each 1/4th segment */
        { 
          iim[L1]=ip[L1]-1+ia[L1]; iip[L1]=iim[L1];
          iim[L2]=ip[L2]-1+ia[L2]; iip[L2]=iim[L2];
          iptm=In4(iim,i4d); iptp=In4(iip,i4d); /* pts either side */
          if (ip[L]==0) { iwp=wherefw[iptp]; iwm=-1; iptm=iptp; }
          else if (ip[L]==i4d[L]) { iwm=wherefw[iptm]; iwp=-1; iptp=iptm; }
          else { iwm=wherefw[iptm]; iwp=wherefw[iptp]; } /* eq no of pts either side */
          if (iwm<0 && iwp<0) continue; /* no equations affected */
          if (iwm==iwp) continue; /* internal surface same equation */
          Loop(i,0,3) id[i]=(ip[i]-1)*2+ia[i];
          geom8areamid(area,xmid,id,L);
          
          if (bij>0) /* reynolds stress terms, use ave stress at points affected */
          { 
            qm=.5*(q[iptm]+q[iptp]);
            if (qm>0)
            {
              rhom=.5*(rho[iptm]+rho[iptp]);
              for (n=0;n<6;n++) 
                uiuj[n2ii[n]][n2jj[n]] = .5*(bij[iptm+iall*n]+bij[iptp+iall*n]); 
              uiuj[0][0]+=1./3.; uiuj[1][1]+=1./3.; uiuj[2][2]+=1./3.;  
              uiuj[1][0]=uiuj[0][1]; uiuj[2][0]=uiuj[0][2]; uiuj[2][1]=uiuj[1][2];
              for (n=0;n<3;n++) for (m=0;m<3;m++) 
                uiuj[n][m]*=2.*qm*qm*rhom; 					
              Loop(n,0,3)   /* each velocity component */
              { 
                drhs=area[0]*uiuj[n][0]+area[1]*uiuj[n][1]+area[2]*uiuj[n][2];
                if (iwm>=0) rhs[n][iwm]-=drhs;
                if (iwp>=0) rhs[n][iwp]+=drhs;
              }
            }
          }
          if (bijm>0 && ip[L]!=0 && ip[L]!=i4d[L]) /* reynolds stress terms, use mid values */
          {  
            mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
            qm=qmid[mmid];
            rhom=rhomid[mmid];
            for (n=0;n<6;n++) 
              uiuj[n2ii[n]][n2jj[n]] = bijm[mmid+iallm*n];
            uiuj[0][0]+=1./3.; uiuj[1][1]+=1./3.; uiuj[2][2]+=1./3.;  
            uiuj[1][0]=uiuj[0][1]; uiuj[2][0]=uiuj[0][2]; uiuj[2][1]=uiuj[1][2];
            for (n=0;n<3;n++) for (m=0;m<3;m++) 
              uiuj[n][m]*=2.*qm*qm*rhom; 					
            Loop(n,0,3)   /* each velocity component */
            {
              drhs=area[0]*uiuj[n][0]+area[1]*uiuj[n][1]+area[2]*uiuj[n][2];
              if (iwm>=0) rhs[n][iwm]-=drhs;
              if (iwp>=0) rhs[n][iwp]+=drhs;
            }
          }
          if (dpdx>0)  /* dpdx */
          { 
            if (iwm>=0)  rhs[0][iwm]-=area[0]*dpdx[0]*xmid[0];
            if (iwp>=0)  rhs[0][iwp]+=area[0]*dpdx[0]*xmid[0];
          }
          if (pg>0) /* pg */
          { 
          }
          if (visc>0 && ip[L]!=0 && ip[L]!=i4d[L]) /* visc but consistent with coefvisc do not do ends */
          { 
            mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));	
            Loop(n,0,3)   /* each momentum component */
            { 
              drhs=area[0]*dudxi[0][mmid+n*iallm]
              +area[1]*dudxi[1][mmid+n*iallm]+area[2]*dudxi[2][mmid+n*iallm];
              drhs*=visc[mmid];
              if (iwm>=0) rhs[n][iwm]+=drhs;
              if (iwp>=0) rhs[n][iwp]-=drhs;
            }
          }
        } /* each 1/4th segment */
		} /* each c.v. */
    } /* each L and ip[3]  */
  }  /* surface integrals */
  /* volume integrals */
  if (zrot!=0)
    Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)
  { 
    if (wherep[In4(ip,i4dp)]<0) continue; /* not valid cont.c.v. */
    geom8volfmid(vol,fmid,ip);
    ig[3]=ip[3];
    Loop3(ia,0,two)  /* corner indices */
    {
      Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
		j=ia[0]+2*ia[1]+4*ia[2];
		ipt[j]=In4(ig,i4d); 
		iw[j]=wherefw[ipt[j]];
    }
    Loop(j,0,8) /* each '1/8' volume */
    {
		if (iw[j]<0) continue;
		Loop(i,0,3) ff[i]=fmid[j][i];
		ffc[0]=(1.-ff[0])*(1.-ff[1])*(1.-ff[2]);
		ffc[1]=ff[0]*(1.-ff[1])*(1.-ff[2]);
		ffc[2]=(1.-ff[0])*ff[1]*(1.-ff[2]);
		ffc[3]=ff[0]*ff[1]*(1.-ff[2]);
		ffc[4]=(1.-ff[0])*(1.-ff[1])*ff[2];
		ffc[5]=ff[0]*(1.-ff[1])*ff[2];
		ffc[6]=(1.-ff[0])*ff[1]*ff[2];
		ffc[7]=ff[0]*ff[1]*ff[2];
		xm[0]=0; xm[1]=0; rum[0]=0; rum[1]=0; rhom=0;
		Loop(i,0,8) 
		{ 
        rhom+=ffc[i]*rho[ipt[j]]; 
        xm[0]+=ffc[i]*x[0][ipt[j]];
        xm[1]+=ffc[i]*x[1][ipt[j]];
        rum[0]+=ffc[i]*rho[ipt[j]]*u[0][ipt[j]];
        rum[1]+=ffc[i]*rho[ipt[j]]*u[1][ipt[j]];
		}
		rhs[0][iw[j]]+= +vol[j]*zrot*2*rum[1];
		rhs[1][iw[j]]+= -vol[j]*zrot*2*rum[0];
		if (nycentrifugal==1)
		{
        rhs[0][iw[j]]+= vol[j]*zrot*rhom*zrot*(xm[0]-xrot);
        rhs[1][iw[j]]+= vol[j]*zrot*rhom*zrot*(xm[1]-yrot);
		}
    }
  }
}	 
