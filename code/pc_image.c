/* contains pc_image */
#include "global.h"
#include "p_plot.h"
/* work with images, combine, output, etc */
/* actions: outjgm, outgif, new, label  */
/* ---------------------------- */
void pc_image(FILE *fpin, FILE *fprint)
{
  char *action, *imagename;
  int *imagetot;
  char *filename,*prename;
  int k=100;
  
  while (k)
  {
    k--;
    action=readname(fpin);
    if (action==NULL) break;
    if(strcmp(action,"end")==0) break;
    if (strcmp(action,"c:")==0) { fseek(fpin,(long)(-3),1);printout("normal","\n");  break; }
    
    if (strcmp(action,"outjgm")==0)
    {
      imagename=readname(fpin);
      filename=readfilename(fpin);
      printout("normal"," writing image %s to file %s, format jgm\n",imagename,filename);
      imagetot=(int*)need(imagename);
      p_outjgm(filename,imagetot+2,imagetot);
      continue;
    }
    if (strcmp(action,"outgif")==0)
    {
      imagename=readname(fpin);
      prename=readfilename(fpin);
      filename=setname2(prename,".gif");     
      printout("normal"," writing image %s to file %s, format gif\n",imagename,filename);
      imagetot=(int*)need(imagename);
      p_outgif(filename,imagetot+2,imagetot);  
      continue;
    }
    if (strcmp(action,"new")==0)
    {
      int ix,iy,i;
      imagename=readname(fpin);
      ix=readint(fpin);
      iy=readint(fpin);
      printout("normal","creating image %s, %d by %d\n",imagename,ix,iy);
      imagetot=(int*)createarray(imagename,2+ix*iy,'i',0);
      imagetot[0]=ix;
      imagetot[1]=iy;
      for (i=2;i<ix*iy+2;i++) imagetot[i]=0;
      continue;
    }
    if (strcmp(action,"label")==0)
    {
      int ifont=1,leftright=0,updown=0,aorp=0,color=1,ix,iy;
      imagename=readname(fpin);
      prename=readname(fpin);
      ix=readint(fpin);
      iy=readint(fpin);
      imagetot=(int*)need(imagename);
      printout("normal","labeling %s with %s centered at %d %d\n",imagename,prename,ix,iy);
      p_symbol(ix,iy,prename,ifont,leftright,updown,aorp,color,imagetot);   
      continue;
    }
    if (strcmp(action,"combine")==0)
    {
      int *imagebase, *imageadd, *image, *imagea; 
      char *imagenamebase;
      int ixll,iyll,i,j,ic; /* add image with lower left corner here in imagebase */
  
      imagenamebase=readname(fpin);
      imagebase=(int*)need(imagenamebase);
      image=imagebase+2;
      imagename=readname(fpin);
      imageadd=(int*)need(imagename);
      imagea=imageadd+2;
      ixll=readint(fpin);
      iyll=readint(fpin);
      printout("normal","add to image %s, image %s, lower-left at %d %d\n",
               imagenamebase,imagename,ixll,iyll);
      for (i=0;i<imageadd[0];i++)
        for (j=0;j<imageadd[1];j++)
        { 
          ic=imagea[i+imageadd[0]*j];
          if (ic==0) continue;
          p_setpixel(i+ixll,j+iyll,ic,image,imagebase);
        }
      continue;
    }
    if (strcmp(action,"recolor")==0)
    {
      int icfrom,icto,i;
      imagename=readname(fpin);
      icfrom=readint(fpin);
      icto=readint(fpin);
      printout("normal","recoloring in %s, from color %d to %d\n",imagename,icfrom,icto);
      imagetot=(int*)need(imagename);
      for (i=2;i<imagetot[0]*imagetot[1]+2;i++)
        if (imagetot[i]==icfrom) imagetot[i]=icto;
    }
    printout("warning pc_image"," action %s not found\n",action);
  }
}
