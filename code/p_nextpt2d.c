/* contains p_nextpt2d */
#include "global.h"
/* cross 2d grid in direction dxy from in and finding out
 solve for four sides solve  xyin + b dyx = xysidei + e (xyside2-xyside1) for b and e
 want e between 0 and 1  b >0  (=0 for collapsed points) */

/* ---------------------------- */
double p_nextpt2d(double dxy[2], double *x[2],int id[2], 
                  int ilowin[2], double fin[2], int ilowout[2], double fout[2])
{
  int di[5]={0,1,1,0,0}, dj[5]={0,0,1,1,0};
  int i,j,m,nc,n;
  double xy[5][2], xyp[2], xypn[2], det,b,e=0;
  double dp,dm;
  int bug=0;
  
  if (bug>0) printout("normal"," in: %d %lg %d %lg\n", ilowin[0],fin[0],ilowin[1],fin[1]);
  /* set 5 points */
  for (n=0;n<5;n++)
  {
    i=di[n]+ilowin[0]; j=dj[n]+ilowin[1];
    m=i+id[0]*j;
    xy[n][0]=x[0][m]; xy[n][1]=x[1][m];
    if (bug>0) printout("normal","n %d xy %lg %lg\n",n,xy[n][0],xy[n][1]);
  }
  /* see if side to omit */
  nc=-1;
  if (fin[0]==0) nc=3;
  else if (fin[0]==1) nc=1;
  else if (fin[1]==0) nc=0;
  else if (fin[1]==1) nc=2;
  
  /* find current point */
  for (i=0;i<2;i++) 
    xyp[i]=(1-fin[0])*(1-fin[1])*x[i][ilowin[0]+ilowin[1]*id[0]]
    +fin[0]*(1-fin[1])*x[i][ilowin[0]+1+ilowin[1]*id[0]]
    +(1-fin[0])*fin[1]*x[i][ilowin[0]+(ilowin[1]+1)*id[0]]
    +fin[0]*fin[1]*x[i][ilowin[0]+1+(ilowin[1]+1)*id[0]];
  if (bug>0) printout("normal","xyp %lg %lg\n",xyp[0],xyp[1]);
  b=-1;
  /* loop over line segments */
  for (n=0;n<4;n++)
  {
    if (n==nc) continue;
    dp=dxy[0]*(xy[n+1][1]-xy[n][1]);
    dm=dxy[1]*(xy[n+1][0]*xy[n][0]);
    det=dxy[0]*(xy[n+1][1]-xy[n][1])-dxy[1]*(xy[n+1][0]-xy[n][0]);
    if (bug>0) printout("normal","n %d det %lg %lg  %lg ",n,dp,dm,det);
    if (det==0) continue;
    dp=dxy[1]*(xy[n][0]-xyp[0])/det;
    dm=dxy[0]*(xy[n][1]-xyp[1])/det;
    e=-(dxy[0]*(xy[n][1]-xyp[1])-dxy[1]*(xy[n][0]-xyp[0]))/det;
    if (bug>0) printout("normal"," e %lg  %lg %lg \n",dp,dm,e);
    if (e<0 || e>1) continue;
    dp=(xy[n][0]-xyp[0])*(xy[n+1][1]-xy[n][1])/det;
    dm=(xy[n][1]-xyp[1])*(xy[n+1][0]-xy[n][0])/det;
    b=((xy[n][0]-xyp[0])*(xy[n+1][1]-xy[n][1])-
       (xy[n][1]-xyp[1])*(xy[n+1][0]-xy[n][0]))/det;
    /* check for roundoff */
    if (abs(b)<1e-8*(abs(dp)+abs(dm))) b=0;
    if (b<0 && nc==-1) continue;
    /* have results */
    if (bug>0) printout("normal","b %lg %lg %lg",dp,dm,b );
    dp=(xy[n][0]-xyp[0]+e*(xy[n+1][0]-xy[n][0]))/dxy[0];
    dm=(xy[n][1]-xyp[1]+e*(xy[n+1][1]-xy[n][1]))/dxy[1];
    if (bug>0) printout("normal"," bx %lg by %lg\n",dp,dm);
    e=max(1.e-6,min(1-1.e-6,e));
    if (n==0) { fout[1]=1; ilowout[1]--; fout[0]=e;}
    else if (n==1) { fout[0]=0; ilowout[0]++; fout[1]=e; }
    else if (n==2) { fout[1]=0; ilowout[1]++; fout[0]=(1-e); }
    else { fout[0]=1; ilowout[0]--; fout[1]=(1-e); }
    break;
  }
  /* check */
  for (i=0;i<2;i++)
  {
  xypn[i]=(1-fout[0])*(1-fout[1])*x[i][ilowout[0]+ilowout[1]*id[0]]
  +fout[0]*(1-fout[1])*x[i][ilowout[0]+1+ilowout[1]*id[0]]
  +(1-fout[0])*fout[1]*x[i][ilowout[0]+(ilowout[1]+1)*id[0]]
  +fout[0]*fout[1]*x[i][ilowout[0]+1+(ilowout[1]+1)*id[0]];
   if (bug>0)  printout("normal","i %d %lg+%lg+%lg+%lg=%lg\n",i,
           (1-fout[0])*(1-fout[1])*x[i][ilowout[0]+ilowout[1]*id[0]],
           fout[0]*(1-fout[1])*x[i][ilowout[0]+1+ilowout[1]*id[0]],
           (1-fout[0])*fout[1]*x[i][ilowout[0]+(ilowout[1]+1)*id[0]],
           fout[0]*fout[1]*x[i][ilowout[0]+1+(ilowout[1]+1)*id[0]],xypn[i]);
  }
  if (bug>0) printout("normal","xpn %lg lhsx %lg rhsx %lg, ypn %lg lhsy %lg rhsy %lg\n",
         xypn[0],xyp[0]+b*dxy[0],xy[n][0]+e*(xy[n+1][0]-xy[n][0]),
         xypn[1],xyp[1]+b*dxy[1],xy[n][1]+e*(xy[n+1][1]-xy[n][1]));
  det=(xypn[1]-xyp[1])*dxy[1]+(xypn[0]-xyp[0])*dxy[0];
  e=(xypn[1]-xyp[1])*dxy[0]-(xypn[0]-xyp[0])*dxy[1];
  if (bug>0) printout("normal","i %d %lg, j %d %lg  dot %lg cross %lg cdd %lg\n",ilowout[0],fout[0],ilowout[1],fout[1],det,e,e/det);
  return b;
}
