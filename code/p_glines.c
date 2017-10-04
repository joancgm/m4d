/*   contains p_glines */
#include "global.h"
#include "p_plot.h"
/* plot grid lines */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

void p_glines(int *gcolor,int linewidth,
              int*i4d, double *x[3],char *clt,
              int irange[3][2],double *gxx[2],double scale,double offset[2],
              int *image, int ipix[2])
{
  int L,L2,L3,ii[4]={0,0,0,0},i,j,ma,mb;
  int gc[3][4],gcmax[3]={0,0,0},color;
  double xp[3],xpa[2],xpb[2];
  
  for (i=0;i<3;i++) for (j=0;j<4;j++) gc[i][j]=gcolor[j+4*i];
  for (i=0;i<3;i++) for (j=0;j<4;j++) gcmax[i]=max(gcmax[i],gc[i][j]);
  for (L=0;L<3;L++)
  { 
    if (gcmax[L]==0) continue;
    if (irange[L][0]==irange[L][1]) continue;
    L2=(L+1)%3; L3=(L+2)%3;
    for (ii[L2]=irange[L2][0];ii[L2]<=irange[L2][1];ii[L2]++)
      for (ii[L3]=irange[L3][0];ii[L3]<=irange[L3][1];ii[L3]++) 
        for (ii[L]=irange[L][0];ii[L]<irange[L][1];ii[L]++)
        {
          ma=In4(ii,i4d);
          ii[L]++;
          mb=In4(ii,i4d);
          ii[L]--;
          if (clt[ma]=='x' || clt[mb]=='x') continue;
          if (clt[ma]=='s' || clt[mb]=='s') color=gc[L][2];
          else if (clt[ma]=='f' || clt[mb]=='f') color=gc[L][0];
          else if (clt[ma]=='w' && clt[mb]=='w') color=gc[L][1];
          else color=gc[L][3];
          if (color==0) continue;
          for (i=0;i<3;i++) xp[i]=x[i][ma];
          for (j=0;j<2;j++) xpa[j]=(Dot(xp,gxx[j])+offset[j])*scale;
          for (i=0;i<3;i++) xp[i]=x[i][mb];
          for (j=0;j<2;j++) xpb[j]=(Dot(xp,gxx[j])+offset[j])*scale;
          p_stlinex(xpa[0],xpb[0],xpa[1],xpb[1],linewidth,color,image,ipix);
        }
  }
}

