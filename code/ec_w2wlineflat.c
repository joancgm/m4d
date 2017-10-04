/* contains c_w2wlineflat */
#include "global.h"
/* calculate the smallest wall-to-wall vector that passes through the midpoint of each continuity control volume,  */
/* set arrays w2wline[noindppts,3], w2wdist[noindppts] */
/* note this only needs to be done once */

/* this starter routine is for flat y,z walls only at j and k grid limits */


#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/*---------------------------c_walldist ------------*/
void ec_w2wlineflat(FILE *fpin, FILE *fprint)
{
  int  *i4d,*whoisp,*noindppts; char *csym;   /* need */
  double *x[3]; char *clt,*cltp;
  double *w2wdist,*w2wline,*w2wnorm; /* create wall distance and wall normal 2 midpoint */
  double *xyzm; /* center location, temporary */
  
  int i,j,m,i4dp[4],iall,ip[4],ig[4];
  double *xc;
  
  double yw[2],zw[2],dy,dz,ys,zs,ydd,zdd;
  
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
  w2wdist=(double *)createarray("w2wdist",noindppts[0],'d',0);
  w2wline=(double *)createarray("w2wline",noindppts[0]*3,'d',0);
  w2wnorm=(double *)createarray("w2wnorm",noindppts[0]*3,'d',0);
  
  yw[0]=x[1][0]; yw[1]=x[1][iall-1];
  zw[0]=x[2][0]; zw[1]=x[2][iall-1];
  printout("normal","yw %lg %lg  zw %lg %lg\n",yw[0],yw[1],zw[0],zw[1]);
  xyzm=(double *)tmalloca(noindppts[0]*3,'d');
  /* initialize wall dist as -1 for real points -2 or -3 for solid or outside pts */
  /* locate the centers for each valid cont c.v. */
  Loop(i,0,noindppts[0]) 
  {
    w2wdist[i]=-1;
    if (cltp[whoisp[i]]=='S' || cltp[whoisp[i]]=='s') w2wdist[i]=-2;
    else
    { 
      iexpand(whoisp[i],i4dp,ip);
      if (ip[0]==0 || ip[1]==0 || ip[2]==0 || 
          ip[0]==i4d[0] || ip[1]==i4d[1] || ip[2]==i4d[2])
      {
        w2wdist[i]=-3; continue;
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
      printout("normal"," ip %d %d %d xyz %lg %lg %lg\n",ip[0],ip[1],ip[2],xc[0],xc[1],xc[2]);
    }
    if (w2wdist[i]<-1.5) continue;
    /* determine wall stuff */
    dy=xc[1]-yw[0]; ys=1;
    if (yw[1]-xc[1] < dy) { dy=yw[1]-xc[1]; ys=-1;}
    dz=xc[2]-zw[0]; zs=1;
    if (zw[1]-xc[2] < dz) { dz=zw[1]-xc[2]; zs=-1;}
    ydd=dz*pow(dy/dz,1./3.);
    zdd=dy*dz/ydd;
    printout("normal","              dy ys %lg %lg dz zs %lg %lg   ydd zdd %lg %lg\n",
           dy,ys,dz,zs,ydd,zdd);
    w2wdist[i]=sqrt((dz+zdd)*(dz+zdd) + (dy+ydd)*(dy+ydd));
    w2wline[i]=0;
    w2wline[i+noindppts[0]]=(dy+ydd)*ys;
    w2wline[i+2*noindppts[0]]=-(dz+zdd)*zs;
    w2wnorm[i]=0;
    w2wnorm[i+noindppts[0]]=(dz+zdd)*ys/w2wdist[i];
    w2wnorm[i+2*noindppts[0]]=(dy+ydd)*zs/w2wdist[i];
    printout("normal","              w2w %lg line %lg %lg %lg norm %lg %lg %lg\n",
           w2wdist[i],w2wline[i],w2wline[i+noindppts[0]],w2wline[i+2*noindppts[0]],
           w2wnorm[i],w2wnorm[i+noindppts[0]],w2wnorm[i+2*noindppts[0]]);
  }
}
