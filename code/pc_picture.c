/* contains pc_picture */
#include "global.h"
#include "p_plot.h"
/* create an 2-d pixel image */
/* input: name, ipix,jpix, pixpergscale, action */

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* ---------------------------- */
void pc_picture(FILE *fpin, FILE *fprint)
{
  char **gnames,**vnames, **fnames;  /* p.items which may be needed */
  double *gxx[2],*grange[3], *vparms, *fparms;
  int *color,*linewidth;
  
  int *i4d; double *x[3], *a[3]; char *clt; /* grid information */
  char *imagename; int ipix[2]; double pixelpergscale; char *action; /* read */
  int *imagetot; /* create */
  int *image; double *gat,*fat; /* alt entries */
  
  int iall,i,j,nystep=0,iadd=0;
  int irange[3][2],iat[4]={0,0,0,0};
  double f[1],offset[2],xat[3];
  double *u[3], *pfill;
  
  /* get needed stuff */
  gnames=(char**)need("p.gnames");
  i4d=(int*)need(gnames[0]);
  iall=Prod4(i4d);
  x[0]=(double*)need(gnames[1]); x[1]=x[0]+iall; x[2]=x[1]+iall;
  a[0]=(double*)need(gnames[2]); a[1]=a[0]+i4d[0]; a[2]=a[1]+i4d[1];
  clt=(char*)need(gnames[3]);
  gxx[0]=(double*)need("p.gxx"); gxx[1]=gxx[0]+3;
  grange[0]=(double*)need("p.grange"); 
  grange[1]=grange[0]+2; grange[2]=grange[1]+2;
  gat=grange[2]+2; fat=gat+3;
  /* read input */
  imagename=readname(fpin);
  ipix[0]=readint(fpin);
  ipix[1]=readint(fpin);
  pixelpergscale=readdouble(fpin);
  action=readname(fpin);
  
  printout("normal","image: %s size %d %d scale %lg action %s\n",
          imagename,ipix[0],ipix[1],pixelpergscale,action);
  /* create image */
  imagetot=(int*)createarray(imagename,2+ipix[0]*ipix[1],'i',0);
  image=imagetot+2;
  imagetot[0]=ipix[0]; imagetot[1]=ipix[1];
  for (i=0;i<ipix[0]*ipix[1];i++) image[i]=0;
  /* set irange  and iat*/
  for (i=0;i<3;i++)
  {
    irange[i][0]=findex(grange[i][0],a[i],i4d[i],f);
    if (f[0]>.99) irange[i][0]++;
    irange[i][1]=findex(grange[i][1],a[i],i4d[i],f);
    if (f[0]>0) irange[i][1]++;
    iat[i]=findex(gat[i],a[i],i4d[i],f);
    if (f[0]>.5) iat[i]++;
    printout("normal","L=%d abc= %lg %lg ijk=%d %d gat %lg iat=%d\n",
            i+1,grange[i][0],grange[i][1],irange[i][0],irange[i][1],gat[i],iat[i]);
  }
  /* set offset */
  for (i=0;i<3;i++) xat[i]=x[i][In4(iat,i4d)];
  /* in general: pixel[j]= ( Dot(xat,gxx[j])+offset[i]) *scale
   located point: pixel[j]= ipix[j]*fat[j] */
  for (j=0;j<2;j++) offset[j]=(fat[j]*ipix[j]/pixelpergscale)-Dot(xat,gxx[j]);
  /* do actions */  
  for (i=0;i<strlen(action);i++)
  {
    if (action[i]=='g') /* grid lines */
    { 
      printout("normal"," draw grid lines %s %c\n",action,action[i]);
      linewidth=(int*)need("p.gparms");
      color=linewidth+1;
      p_glines(color,linewidth[0],i4d,x,clt,irange,gxx,pixelpergscale,offset,image,ipix);
    }
    else if (action[i]=='v') /* velocity vectors */
    {
      printout("normal"," draw velocity vectors %s %c\n",action,action[i]);
      vnames=(char **)need("p.vnames");
      vparms=(double*)need("p.vparms");
      for (j=0;j<3;j++) u[j]=(double*)need(vnames[j]);
      p_vectors(u,vparms,i4d,x,clt,irange,gxx,pixelpergscale,offset,image,ipix);
    }
    else if (action[i]=='f') /* contour fill */
    {
      fnames=(char **)need("p.fnames");
      printout("normal"," fill contour %s %c\n",action,action[i]);
      fparms=(double*)need("p.fparms");
      pfill=(double*)need(fnames[0]);
      j=arraysize(fnames[0]);
      
      if (j==iall) { nystep=0; iadd=0;}
      else if (j==(i4d[0]-1)*(i4d[1]-1)*(i4d[2]-1)*i4d[3]) { nystep=1; iadd=0;}
      else if (j==(i4d[0]+1)*(i4d[1]+1)*(i4d[2]+1)*i4d[3]) { nystep=1; iadd=1;}
      else
      {
        printout("error pc_picture","incompatible size, not grid, g+1 or g-1, cannot contour %s, exit\n",
                fnames[0]);
        exitm4d(0);
      }
      printout("normal"," fill with prop %s, size g+%d method %d (0=linear, 1=stepwise)\n",
              fnames[0],iadd*2-nystep,nystep);
      p_fillcontour(pfill,nystep,iadd,fnames,fparms,
                    i4d,x,clt,irange,gxx,pixelpergscale,offset,image,ipix);
    }
  }
}
