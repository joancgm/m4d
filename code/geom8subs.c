/* contains geom8init, geom8vol, geom8volfmid, geom8fx27, geom8areamid, geom8gradvol, geom8gradvxf, c_geom8print */
#include "global.h"
/* routines to calculate "1/8th" volumes or "1/4th" areas  based on x and cvd  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* for cvdcdouble */
#define Mi(i,j,k,n) ((j)+2*i4d[1]*((k)+2*i4d[2]*((i)+i4d[0]*n)))
#define Mj(i,j,k,n) ((k)+2*i4d[2]*((i)+2*i4d[0]*((j)+i4d[1]*n)))
#define Mk(i,j,k,n) ((i)+2*i4d[0]*((j)+2*i4d[1]*((k)+i4d[2]*n)))
/* statics because routines called for individual points */
static int *i4d,i4dd[4],iall,ialld;
static double *x[3],*xd[3],*cvdc,*cvdci,*cvdcj,*cvdck;
/* --------------------------------*/
void geom8init(void)  /* locate arrays  for later use */
{ 
  int i;
  i4d=(int *)need("idim4d");   
  iall=Prod4(i4d);
  Loop(i,0,3) i4dd[i]=2*i4d[i]-1; i4dd[3]=i4d[3];
  ialld=Prod4(i4dd);
  x[0]=(double *)need("xyz"); 
  x[1]=x[0]+iall; 
  x[2]=x[1]+iall;
  cvdc=(double *)need("cvdc");
  cvdci=(double *)need("cvdcdouble");
  cvdcj=cvdci+4*iall; cvdck=cvdcj+4*iall;
  xd[0]=(double *)find("xyzdouble");
  xd[1]=xd[0]+ialld;
  xd[2]=xd[1]+ialld;
}
/*--------------------------------------*/
/* set   vol = the 8 1/8th volumes for cont c.v. at ip */
void geom8vol(double *vol, int *ip)
{ 
  int i,j,k,ia[3],ib[3];
  double f[27][3],xg[27][3],dx[3][3],da[3];
  int two[3]={2,2,2};
  
  geom8fx27(f,xg,ip);
  Loop3(ia,0,two)
  {
    Loop(i,0,3) Loop(j,0,3) dx[i][j]=0;
    Loop3(ib,0,two)
    { 
      j=ia[0]+ib[0]+3*(ia[1]+ib[1])+9*(ia[2]+ib[2]);
      Loop(i,0,3) Loop(k,0,3) dx[i][k]+=.25*xg[j][k]*(2*ib[i]-1);
    }
    Cross(dx[0],dx[1],da);
    j=ia[0]+2*ia[1]+4*ia[2];
    vol[j]=Dot(dx[2],da);
    /*	  Loop(i,0,3)   printout("normal","dxL %d %lg %lg %lg ",i,dx[i][0],dx[i][1],dx[i][2]);
     printout("normal","da %lg %lg %lg   vol  %d %lg\n",da[0],da[1],da[2],j,vol[j]); */
  }
}
/*--------------------------------------*/
/* set   vol = the 8 1/8th volumes for cont c.v. at ip, and fmid the interpolation for the center*/
void geom8volfmid(double *vol, double (*fmid)[3], int *ip)
{
  int i,j,k,n,ia[3],ib[3];
  double f[27][3],xg[27][3],dx[3][3],da[3];
  int two[3]={2,2,2};
  
  geom8fx27(f,xg,ip);
  Loop3(ia,0,two)
  {
    Loop(i,0,3) Loop(j,0,3) dx[i][j]=0;
    n=ia[0]+2*ia[1]+4*ia[2];
    Loop(j,0,3) fmid[n][j]=0;
    Loop3(ib,0,two)
    {
      j=ia[0]+ib[0]+3*(ia[1]+ib[1])+9*(ia[2]+ib[2]);
      Loop(i,0,3) Loop(k,0,3) dx[i][k]+=.25*xg[j][k]*(2*ib[i]-1);
      Loop(k,0,3) fmid[n][k]+=.125*f[j][k];
    }
    Cross(dx[0],dx[1],da);
    vol[n]=Dot(dx[2],da);
  }
  /*	  printout("normal","geom8volfmid  ip %d%d %d \n",ip[0],ip[1],ip[2]);
   Loop(i,0,8) printout("normal","vol f %lg %lg %lg %lg\n",vol[i],fmid[i][0],fmid[i][1],fmid[i][2]); */
}

/* set  the 3x3x3 interpolation factors and the 3x3x3 grid for the subdivided cont c.v. */
/* redone for cvdc and cvdc3 */
void geom8fx27(double (*f)[3], double (*xg)[3], int *ip)
{ 
  int i,j,ia[3],ib[3],L,L2,L3,iglow[4],ig[4],ipt,iptc[8],ja,jb;
  int two[3]={2,2,2};
  int three[3]={3,3,3};
  int ifac[3]={1,3,9};
  
  Loop3(ib,0,three)   /* initial uniform f */
  {
    j=ib[0]+3*ib[1]+9*ib[2];
	 Loop(i,0,3) f[j][i]=.5*ib[i];
  }
  /* xyz at 3x3x3 grid corners  and indices for them */
  Loop(i,0,3) iglow[i]=ip[i]-1; iglow[3]=ip[3];  
  ig[3]=ip[3];
  Loop3(ia,0,two)   
  { 
    Loop(i,0,3) ig[i]=iglow[i]+ia[i];
	 ipt=In4(ig,i4d);
	 iptc[ia[0]+2*ia[1]+4*ia[2]]=ipt;
	 j=2*ia[0]+3*2*ia[1]+9*2*ia[2];
	 Loop(i,0,3) xg[j][i]=x[i][ipt];  
  }
  /* f then x  based on cvdcdouble */
  /* complete f */
  Loop(ja,0,3) Loop(jb,0,3)
  { 
    j=1+3*ja+9*jb;
    f[j][0]=cvdci[Mi(iglow[0],2*iglow[1]+ja,2*iglow[2]+jb,iglow[3])];
    j=ja+3+9*jb;
    f[j][1]=cvdcj[Mj(2*iglow[0]+ja,iglow[1],2*iglow[2]+jb,iglow[3])];
    j=ja+3*jb+9;
    f[j][2]=cvdck[Mk(2*iglow[0]+ja,2*iglow[1]+jb,iglow[2],iglow[3])];
  }
  Loop(L,0,3)    /* lines */
  {
    L2=(L+1)%3; 
	 L3=(L+2)%3;
	 Loop(ia[L2],0,2) Loop(ia[L3],0,2) 
    {
      ia[L]=0;
      Loop(i,0,3) ib[i]=2*ia[i]; 
      ib[L]=1; 
      j=ib[0]+3*ib[1]+9*ib[2];
      Loop(i,0,3) xg[j][i]=(1.-f[j][L])*xg[j-ifac[L]][i] +f[j][L]*xg[j+ifac[L]][i];
    }
  }
  Loop(L,0,3)     /* sides */
  {
    L2=(L+1)%3;  L3=(L+2)%3;
	 for (ib[L]=0;ib[L]<3;ib[L]+=2)
    { 
      ib[L2]=1; ib[L3]=1; 
      j=ib[0]+3*ib[1]+9*ib[2];
      Loop(i,0,3) xg[j][i]=(1.-f[j][L2])*(1-f[j][L3])*xg[j-ifac[L2]-ifac[L3]][i]
      +(1.-f[j][L2])*(f[j][L3])*xg[j-ifac[L2]+ifac[L3]][i]
      +(f[j][L2])*(1-f[j][L3])*xg[j+ifac[L2]-ifac[L3]][i]
      +(f[j][L2])*(f[j][L3])*xg[j+ifac[L2]+ifac[L3]][i];
    }
  }
  j=13;    /* center */
  Loop(i,0,3)   /* center */
  xg[j][i]=(1-f[j][0])*(1-f[j][1])*(1-f[j][2])*xg[j-ifac[0]-ifac[1]-ifac[2]][i]
  +(f[j][0])*(1-f[j][1])*(1-f[j][2])*xg[j+ifac[0]-ifac[1]-ifac[2]][i]
  +(1-f[j][0])*(f[j][1])*(1-f[j][2])*xg[j-ifac[0]+ifac[1]-ifac[2]][i]
  +(1-f[j][0])*(1-f[j][1])*(f[j][2])*xg[j-ifac[0]-ifac[1]+ifac[2]][i]
  +(f[j][0])*(f[j][1])*(1-f[j][2])*xg[j+ifac[0]+ifac[1]-ifac[2]][i]
  +(1-f[j][0])*(f[j][1])*(f[j][2])*xg[j-ifac[0]+ifac[1]+ifac[2]][i]
  +(f[j][0])*(1-f[j][1])*(f[j][2])*xg[j+ifac[0]-ifac[1]+ifac[2]][i]
  +(f[j][0])*(f[j][1])*(f[j][2])*xg[j+ifac[0]+ifac[1]+ifac[2]][i];
  
}
/*-------------------------------------*/
void geom8areamid(double *area, double *xmid, int *id, int L)
{ 
  int i,j,L2,L3,ia[3],ii[4],ib,ic,m,is[3];
  double dx[2][3];
  L2=(L+1)%3; L3=(L+2)%3;
  Loop(i,0,3) {area[i]=0; xmid[i]=0;}
  if (xd[0]==0) 
  {
    printout("error geom8areamid","Error: need xyzdouble for geom8areamid\n"); 
    exitm4d(0);
  }
  Loop(i,0,4) {if (id[i]<0) return; }  /* check out of bounds */
  if (id[L]>i4dd[L]-1) return;
  if (id[L2]>i4dd[L2]-2) return;
  if (id[L3]>i4dd[L3]-2) return;
  if (id[3]>i4dd[3]-1) return;
  
  ii[3]=id[3];  ia[L]=0; is[L]=1;
  Loop(i,0,2) Loop(j,0,3) dx[i][j]=0;
  
  Loop(ib,0,2) Loop(ic,0,2)
  { 
    ia[L2]=ib; ia[L3]=ic;
	 is[L2]=2*ib-1; is[L3]=2*ic-1;
	 Loop(i,0,3) ii[i]=id[i]+ia[i];
	 m=In4(ii,i4dd);
	 Loop(i,0,3) 
    {
      dx[0][i]+=.5*xd[i][m]*is[L2];
      dx[1][i]+=.5*xd[i][m]*is[L3];
      xmid[i]+=.25*xd[i][m];
    }
  }
  /*	 printout("normal","dx0 %lg %lg %lg dx1 %lg %lg %lg\n",
   dx[0][0],dx[0][1],dx[0][2],dx[1][0],dx[1][1],dx[1][2]); */
  Cross(dx[0],dx[1],area);
}
/*----------------------------*/
/* for the 8  "1/8" volumes  calc grad dVol by an area integral  given cont c.v. index */
void geom8gradvol(double (*gradv)[8][3], int *ip)
{ 
  double f[27][3],x[27][3];
  geom8fx27(f,x,ip);
  geom8gradvxf(gradv,x,f);
}

void geom8gradvxf(double (*gradv)[8][3], double (*x)[3], double (*f)[3])
{
  int i,j,k,n,ia[3],ib[3],ic[3],L,L2,L3;
  double fmid[3][2], dx[3][3],da[3];
  int two[3]={2,2,2};
  
  Loop(i,0,8) Loop(j,0,8) Loop(k,0,3) gradv[i][j][k]=0;
  
  Loop3(ia,0,two)  /* each 1/8th */
  {
    j=ia[0]+3*ia[1]+9*ia[2];
	 n=ia[0]+2*ia[1]+4*ia[2];
    Loop(L,0,3)
    {
      L2=(L+1)%3; L3=(L+2)%3;
      Loop(ib[L],0,2)
      { 
        Loop(i,0,3) fmid[i][1]=0;
        Loop(i,0,3) Loop(k,0,3) dx[i][k]=0;
        Loop(ib[L2],0,2) Loop(ib[L3],0,2)
        {
          k=j+ib[0]+3*ib[1]+9*ib[2];
          Loop(i,0,3)
          {
            fmid[i][1]+=.25*f[k][i];
            dx[L2][i]+=.5*x[k][i]*(2*ib[L2]-1);
            dx[L3][i]+=.5*x[k][i]*(2*ib[L3]-1);
          }
        }
        /*		 printout("normal","L %d  ib[L] %d dxL2 %lg %lg %lg dxL3 %lg %lg %lg fmid %lg %lg %lg\n",
         L,ib[L],dx[L2][0],dx[L2][1],dx[L2][2],dx[L3][0],dx[L3][1],dx[L3][2], 
         fmid[0][1],fmid[1][1],fmid[2][1]); */
        Cross(dx[L2],dx[L3],da);
        Loop(i,0,3) 
        {
          da[i]=da[i]*(2*ib[L]-1); 
          fmid[i][0]=1.-fmid[i][1]; 
        }
        Loop3(ic,0,two)
        {
          k=ic[0]+2*ic[1]+4*ic[2];
          Loop(i,0,3) gradv[n][k][i]+=da[i]*fmid[0][ic[0]]*fmid[1][ic[1]]*fmid[2][ic[2]];
        }
      }
    }
  }
}
/*--------------------------------------*/ 
/* a command routine to do print tests */
void c_geom8print(FILE *fpin, FILE *fprint)
{ 
  int i,j,n,id[4],is[4],ie[4],ips[4],ipe[4];
  double area[3],xmid[3],vol[8],f[27][3],x[27][3],gradv[8][8][3];
  
  geom8init();
  Loop(i,0,4) 
  {
    is[i]=readint(fpin); is[i]=max(0,min(i4dd[i]-1,is[i]));
	 ie[i]=readint(fpin)+1; ie[i]=max(1,min(i4dd[i],ie[i]));
  }
  if (xd[0]>0)
  {
    Loop(id[3],is[3],ie[3]) Loop(id[2],is[2],ie[2])
    Loop(id[1],is[1],ie[1]) Loop(id[0],is[0],ie[0]) 
    { 
      printout("print","g8 id %d %d %d %d \n",id[0],id[1],id[2],id[3]);
      Loop(i,0,3)  
      {
        geom8areamid(area,xmid,id,i); 
        printout("print","   L %d area %lg %lg %lg xmid %lg %lg %lg\n", 
                i,area[0],area[1],area[2],xmid[0],xmid[1],xmid[2]); 
      } 
    } 
  }
  else printout("print"," xd not available for geom8areamid\n"); 
  
  Loop(i,0,3) 
  {
    ips[i]=max(1,min(i4d[i],is[i]/2));
    ipe[i]=max(2,min(i4d[i]+1,ie[i]/2));	
  }
  Loop(id[3],is[3],ie[3]) Loop(id[2],ips[2],ipe[2])
  Loop(id[1],ips[1],ipe[1]) Loop(id[0],ips[0],ipe[0])
  {
    printout("print","g8 ip= %d %d %d %d\n",id[0],id[1],id[2],id[3]);	
	 geom8fx27(f,x,id);
    Loop(n,0,3) 
    { 
      printout("print","f %d",n);
      Loop(i,0,27)
      {printout("print"," %lg",f[i][n]); if ((i+1)%9==0) printout("print","\n"); }
      printout("print","x %d",n);
      Loop(i,0,27)
      {printout("print"," %lg",x[i][n]); if ((i+1)%9==0) printout("print","\n"); }
    } 
    geom8vol(vol,id);
	 printout("print","vol");
	 Loop(i,0,8) printout("print"," %lg",vol[i]); printout("print","\n");
	 
	 geom8gradvol(gradv,id);
	 Loop(i,0,8) Loop(n,0,3)
    {
      printout("print","qtr %d grad %d  ",i,n);
      Loop(j,0,8) printout("print"," %lg",gradv[i][j][n]);
      printout("print","\n");
    } 
  }
}
