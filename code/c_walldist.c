/* contains c_walldist */
#include "global.h"
/* calculate distance to the neareat wall for the midpoint of each continuity control volume */
/* set array walldist[noindppts], note this only needs to be done once */


#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/*---------------------------c_walldist ------------*/
void c_walldist(FILE *fpin, FILE *fprint)
{
  int  *i4d,*whoisp,*noindppts; char *csym;   /* need */
  double *x[3]; char *clt,*cltp;
  double *walldist,*walln2m; /* create wall distance and wall normal 2 midpoint */
  double *xyzm; /* center location, temporary */
  
  int i,j,n,m,i4dp[4],iall,ip[4],ig[4],ipw,ia[3],m8[8],m4[5],ntri;
  int two[3]={2,2,2};
  int nrep=0,ib,na,nt;
  double *xc,aa[3][3],vnorm,dxrep[3][3],xp[3],xt[3],xxp,xxt,f[4],bb[3][3],wd,det;
  double v1t[3][3],v2t[3][3],vn[3],*vnnew;
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  whoisp=(int *)need("whoisp");
  noindppts=(int *)need("noindppts");
  csym=(char *)need("csym");
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[1]+iall;
  clt=(char *)need("clt");
  cltp=(char*)need("cltp");
  walldist=(double *)createarray("walldist",noindppts[0],'d',0);
  walln2m=(double *)createarray("walln2m",noindppts[0]*3,'d',0);
  
  xyzm=(double *)tmalloca(noindppts[0]*3,'d');
  /* initialize wall dist as -1 for real points -2 or -3 for solid or outside pts */
  /* locate the centers for each valid cont c.v. */
  Loop(i,0,noindppts[0]) 
  {
    walldist[i]=-1;
    if (cltp[whoisp[i]]=='S' || cltp[whoisp[i]]=='s') walldist[i]=-2;
    else
    { 
      iexpand(whoisp[i],i4dp,ip);
      if (ip[0]==0 || ip[1]==0 || ip[2]==0 || 
          ip[0]==i4d[0] || ip[1]==i4d[1] || ip[2]==i4d[2])
      {
        walldist[i]=-3; continue;
      }
      xc=xyzm+i*3;
      ig[3]=ip[3];
      Loop(j,0,3) xc[j]=0;
      for (ig[0]=ip[0]-1;ig[0]<ip[0]+1;ig[0]++)
        for (ig[1]=ip[1]-1;ig[1]<ip[1]+1;ig[1]++)
          for (ig[2]=ip[2]-1;ig[2]<ip[2]+1;ig[2]++)
          { 
            m=In4(ig,i4d);
            Loop(j,0,3) xc[j]+=.125*x[j][m];
          }	
    }
  }
  /* see if need to check repeats */
  if (csym[6]=='r') {Loop(j,0,3) dxrep[0][j]=x[j][i4d[0]-1]-x[j][0]; nrep++; }
  if (csym[8]=='r') {Loop(j,0,3) dxrep[nrep][j]=x[j][i4d[0]*(i4d[1]-1)]-x[j][0]; nrep++; }
  if (csym[10]=='r') {Loop(j,0,3) dxrep[nrep][j]=x[j][i4d[0]*i4d[1]*(i4d[2]-1)]-x[j][0]; nrep++; }
  
  /* look for wall segments, on points about p-points type F */
  Loop(ipw,0,noindppts[0])
  {  
    /* iexpand(whoisp[ipw],i4dp,ip); */
    if (cltp[whoisp[ipw]]!='F') continue;
    if (walldist[ipw]<-1.5) continue;
    iexpand(whoisp[ipw],i4dp,ip);
    /* printout("normal","checking wall for ip %d %d %d %d",ip[0],ip[1],ip[2],ip[3]); */
    Loop3(ia,0,two)
    { 
      Loop(j,0,3) ig[j]=ip[j]-1+ia[j];
      j=ia[0]+2*ia[1]+4*ia[2];
      m8[j]=In4(ig,i4d);
      /* printout("normal"," m %d clt %c",m8[j],clt[m8[j]]); */
    }
    /* printout("normal","\n"); */
    Loop(ib,0,6)  /* set up 4 ordered possible wall points */
    { 
      if (ib==0) {m4[0]=m8[0]; m4[1]=m8[2]; m4[2]=m8[6]; m4[3]=m8[4]; }
      else if (ib==1) {m4[0]=m8[1]; m4[1]=m8[5]; m4[2]=m8[7]; m4[3]=m8[3]; }
      else if (ib==2) {m4[0]=m8[0]; m4[1]=m8[4]; m4[2]=m8[5]; m4[3]=m8[1]; }
      else if (ib==3) {m4[0]=m8[2]; m4[1]=m8[3]; m4[2]=m8[7]; m4[3]=m8[6]; }
      else if (ib==4) {m4[0]=m8[0]; m4[1]=m8[1]; m4[2]=m8[3]; m4[3]=m8[2]; }
      else if (ib==5) {m4[0]=m8[4]; m4[1]=m8[6]; m4[2]=m8[7]; m4[3]=m8[5]; }
      if (clt[m4[0]]!='w' || clt[m4[1]]!='w' || clt[m4[2]]!='w' || clt[m4[3]]!='w') continue;
      m4[4]=m4[0];  /* so can check 2 triangles starting 0,1,2  then 2,3,4 */
      for (ntri=0; ntri<3; ntri+=2) /* each start for 3 pt triangle */
      { 
        /* make all relative to first point, set up equations f[1]*v1+f[2]*v2+f[3]*vn=p  */
        Loop(n,0,3) 
        {
          aa[0][n]=x[n][m4[ntri+1]]-x[n][m4[ntri]];
          aa[1][n]=x[n][m4[ntri+2]]-x[n][m4[ntri]];
        }
        if (Sqr(aa[0])==0 || Sqr(aa[1])==0) continue;
        Cross(aa[0],aa[1],aa[2]);
        vnorm=sqrt(Sqr(aa[2]));
        Loop(n,0,3) aa[2][n]/=vnorm;  /* normal vector for triangle */
        det=aa[0][0]*(aa[1][1]*aa[2][2]-aa[1][2]*aa[2][1])
        +aa[1][0]*(aa[2][1]*aa[0][2]-aa[2][2]*aa[0][1])
        +aa[2][0]*(aa[0][1]*aa[1][2]-aa[0][2]*aa[1][1]);
        
        /* loop over each midpoint, rather inefficient but complete */
        Loop(i,0,noindppts[0]) if (walldist[i]>=-1.5)
        { 
          /* iexpand(whoisp[i],i4dp,ipm); */
          xc=xyzm+i*3;
          Loop(n,0,3) xp[n]=xc[n]-x[n][m4[ntri]];
          if (nrep>0) /* check if closest to init point is across repeat boundary */
          {
            xxp=Sqr(xp);
            Loop(j,0,nrep)
            { 
              Loop(n,0,3) xt[n]=xp[n]+dxrep[j][n];
              xxt=Sqr(xt);
              if (xxt<xxp) {xxp=xxt; xp[0]=xt[0]; xp[1]=xt[1]; xp[2]=xt[2]; continue; }
              Loop(n,0,3) xt[n]=xp[n]-dxrep[j][n];
              xxt=Sqr(xt);
              if (xxt<xxp) {xxp=xxt; xp[0]=xt[0]; xp[1]=xt[1]; xp[2]=xt[2];}
            }
          } /* if nrep */
          for (nt=0;nt<3;nt++) for (na=0;na<3;na++)  bb[na][nt]=aa[na][nt];
          for (nt=0;nt<3;nt++) bb[2][nt]=xp[nt];
          f[2]=(bb[0][0]*(bb[1][1]*bb[2][2]-bb[1][2]*bb[2][1])
                +bb[1][0]*(bb[2][1]*bb[0][2]-bb[2][2]*bb[0][1])
                +bb[2][0]*(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1]))/det;
          if (f[2]<0) continue; /* wrong side of surface */
          for (nt=0;nt<3;nt++) {bb[2][nt]=aa[2][nt]; bb[0][nt]=xp[nt]; }
          f[0]=(bb[0][0]*(bb[1][1]*bb[2][2]-bb[1][2]*bb[2][1])
                +bb[1][0]*(bb[2][1]*bb[0][2]-bb[2][2]*bb[0][1])
                +bb[2][0]*(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1]))/det;
          for (nt=0;nt<3;nt++) {bb[0][nt]=aa[0][nt]; bb[1][nt]=xp[nt]; }
          f[1]=(bb[0][0]*(bb[1][1]*bb[2][2]-bb[1][2]*bb[2][1])
                +bb[1][0]*(bb[2][1]*bb[0][2]-bb[2][2]*bb[0][1])
                +bb[2][0]*(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1]))/det;
          if (f[0]>0 && f[1]>0 && 1-f[0]-f[1]>0) { j=-1; wd=f[2]; vnnew=aa[2]; }
          else 	/* look for closest normal point on edge of triangle */
          { 
            xxp=Sqr(xp); vnnew=xp;
            for (n=0;n<3;n++) 
            {
              v1t[0][n]=aa[0][n]; v2t[0][n]=xp[n];
              v1t[1][n]=aa[1][n]-aa[0][n]; v2t[1][n]=xp[n]-aa[0][n];
              v1t[2][n]=aa[1][n]; v2t[2][n]=xp[n];
            }
            j=0; for (na=0;na<3;na++) 
            { 
              f[na]=Dot(v1t[na],v2t[na])/Dot(v1t[na],v1t[na]);
              if (f[na]>=0 && f[na]<=1.)
              {  
                for (n=0;n<3;n++) vn[n]=v2t[na][n]-f[na]*v1t[na][n];
                if (Sqr(vn)<xxp) { j=na+1; xxp=Sqr(vn); vnnew=vn; }
              }
            }
            if (j==0) /* check other points */
            {
              if (Sqr(v2t[1])<xxp) {xxp=Sqr(v2t[1]); vnnew=v2t[1]; }
              for (n=0;n<3;n++) v2t[2][n]=xp[n]-aa[1][n];
              if (Sqr(v2t[2])<xxp) {xxp=Sqr(v2t[2]); vnnew=v2t[2]; }
            }
            wd=sqrt(xxp);
          }
          if (walldist[i]<0 || wd<walldist[i])
          {  
            walldist[i]=wd;
            vnorm=sqrt(Sqr(vnnew));
            if (vnorm>0) Loop(n,0,3) vnnew[n]/=vnorm;
            Loop(n,0,3) walln2m[i+n*noindppts[0]]=vnnew[n];
          }
        } /* loop over mid-points */
      } /* each triangle */
    } /* ib ordered wall points */
  } /* eack wall segment */
  printout("normal"," set walldist[noindppts], distance from mid-c.v. to nearest wall\n");
}	
