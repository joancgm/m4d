/*   contains p_vectors */
#include "global.h"
#include "p_plot.h"
/* plot vectors */
/* vparms[11]: dt head headf nybet lwpixel colotz colors colore dUz Uns Une */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Vabs(a) (sqrt((double)( *(a) * *(a) + *((a)+1) * *((a)+1) + *((a)+2) * *((a)+2) )))

void p_vectors(double *u[3], double *vparms,
              int*i4d, double *x[3],char *clt,
              int irange[3][2],double *gxx[2],double scale,double offset[2],
              int *image, int ipix[2])
{
  double xbase[3],xtip[3],xnorm[3],upt[3],dd;
  int i,j,k,ii[4]={0,0,0,0},m[8],nave;
  /* int ipa[2],ipb[2]; */
  int nybet[3],linewidth,color,nynorm;
  double dt,headmax,headf,dUz,Uns,Une,colorz,colors,colore,dn;
  double xplbase[2],xpltip[2],head,headx,heady,pa[2],pb[2];
  
  dt=vparms[0];
  headmax=vparms[1];
  headf=vparms[2];
  for (i=0;i<3;i++)
  {
    nybet[i]=0;
    if (vparms[3] !=0 && irange[i][1]>irange[i][0]) nybet[i]=1;
  }
  linewidth=vparms[4]+.01;
  colorz=vparms[5];
  colors=vparms[6];
  colore=vparms[7];
  dUz=vparms[8];
  Uns=vparms[9];
  Une=vparms[10];
  nynorm=0; if (colore != colors && Uns != Une) nynorm=1;
  /* printout("normal","range %d %d %d %d %d %d\n",
         irange[0][0],max(irange[0][0],irange[0][1]-nybet[0]),
         irange[1][0],max(irange[1][0],irange[1][1]-nybet[1]),
         irange[2][0],max(irange[2][0],irange[2][1]-nybet[2])); */
  { 
    for (ii[0]=irange[0][0];ii[0]<=max(irange[0][0],irange[0][1]-nybet[0]);ii[0]++)
      for (ii[1]=irange[1][0];ii[1]<=max(irange[1][0],irange[1][1]-nybet[1]);ii[1]++)
        for (ii[2]=irange[2][0];ii[2]<=max(irange[2][0],irange[2][1]-nybet[2]);ii[2]++)
        {
          /* set list of points to be averaged */
          nave=0;
          for (i=0;i<2;i++)
          {
            if (i>0 && nybet[0]==0) continue;
            for (j=0;j<2;j++) 
            {
              if (j>0 && nybet[1]==0) continue;
              for (k=0;k<2;k++)
              {
                if (k>0 && nybet[2]==0) continue;
                ii[0]+=i; ii[1]+=j; ii[2]+=k;
                m[nave]=In4(ii,i4d);
                ii[0]-=i; ii[1]-=j; ii[2]-=k;
                nave++;
              }
            }
          }
          /* average points and vectors */
          dd=nave;
          dd=1/dd;
          for (i=0;i<3;i++) { xbase[i]=0; upt[i]=0; }
          for (i=0;i<nave;i++)
          {
            for (j=0;j<3;j++) 
            {
              xbase[j]+=dd*x[j][m[i]];
              upt[j]+=dd*u[j][m[i]];
            }
          } 
        /* printout("normal","ii %d %d %d nave %d upt %lg %lg %lg\n",
                 ii[0],ii[1],ii[2],nave,upt[0],upt[1],upt[2]); */
          if (upt[0]==0 && upt[1]==0 && upt[2]==0) continue;
          for (i=0;i<3;i++) xtip[i]=xbase[i]+dt*upt[i];
          for (i=0;i<2;i++) 
          {
            xplbase[i]=Dot(xbase,gxx[i]);
            xpltip[i]=Dot(xtip,gxx[i]);              
          }
          dd=sqrt((xpltip[1]-xplbase[1])*(xpltip[1]-xplbase[1])+
                  (xpltip[0]-xplbase[0])*(xpltip[0]-xplbase[0]));
          if (dd==0) continue;
          head=min(dd*headf,headmax);
          headx=(xpltip[0]-xplbase[0])*head/dd;
          heady=(xpltip[1]-xplbase[1])*head/dd;
         /* printout("normal","dd,head %lg %lg %lg %lg\n",
                 dd,head,headx,heady); */
          /* determine color */
          if (nynorm==0) color=colorz+.01;
          else /* determine normal to get color */
          {
            Cross(gxx[1],gxx[0],xnorm);
            dn=Vabs(xnorm);
            dd=Dot(xnorm,upt)/dn;
            if (abs(dd)<dUz) color=colorz+.01;
            else if (dd<=Uns) color=colors+.01;
            else if (dd>=Une) color=colore+.01;
            else color=colors+.01+(colore-colors)*(dd-Uns)/(Une-Uns);
           /* printout("normal","xnorm %lg %lg %lg dd %lg U3 %lg color %d \n",
                   xnorm[0]/dn,xnorm[1]/dn,xnorm[2]/dn,dd,upt[2],color); */
          }
          for (j=0;j<2;j++) pa[j]=(Dot(xtip,gxx[j])+offset[j])*scale;
          for (j=0;j<2;j++) pb[j]=(Dot(xbase,gxx[j])+offset[j])*scale;
          /*printout("normal","pa pb %lg %lg %lg %lg \n",pa[0],pa[1],pb[0],pb[1]); */
          p_stlinex(pa[0],pb[0],pa[1],pb[1],linewidth,color,image,ipix);
          pb[0]=pa[0]+(-headx-.4*heady)*scale;
          pb[1]=pa[1]+(-heady+.4*headx)*scale;
          p_stlinex(pa[0],pb[0],pa[1],pb[1],linewidth,color,image,ipix);
          pb[0]=pa[0]+(-headx+.4*heady)*scale;
          pb[1]=pa[1]+(-heady-.4*headx)*scale;
          p_stlinex(pa[0],pb[0],pa[1],pb[1],linewidth,color,image,ipix);
          
          
         /* for (j=0;j<2;j++) ipa[j]=(Dot(xtip,gxx[j])+offset[j])*scale;
          for (j=0;j<2;j++) ipb[j]=(Dot(xbase,gxx[j])+offset[j])*scale;
          printout("normal","ipa ipb %d %d %d %d ",ipa[0],ipa[1],ipb[0],ipb[1]);
          p_stline(ipa[0],ipb[0],ipa[1],ipb[1],linewidth,color,image,ipix);
          ipb[0]=ipa[0]+(-headx-.4*heady)*scale;
          ipb[1]=ipa[1]+(-heady+.4*headx)*scale;
          printout("normal","ipa ipb %d %d %d %d ",ipa[0],ipa[1],ipb[0],ipb[1]);
          p_stline(ipa[0],ipb[0],ipa[1],ipb[1],linewidth,color,image,ipix);
          ipb[0]=ipa[0]+(-headx+.4*heady)*scale;
          ipb[1]=ipa[1]+(-heady-.4*headx)*scale;
          printout("normal","ipa ipb %d %d %d %d\n ",ipa[0],ipa[1],ipb[0],ipb[1]);
          p_stline(ipa[0],ipb[0],ipa[1],ipb[1],linewidth,color,image,ipix);
          */
        }
  }
}


