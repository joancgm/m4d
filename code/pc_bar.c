/* contains pc_bar */
#include "global.h"
#include "p_plot.h"
/* create a color bar image */
/* input: name, ipix,jpix, label */

/* ---------------------------- */
void pc_bar(FILE *fpin, FILE *fprint)
{
  char **fnames;  /* p.items needed */
  double *fparms;
  int *fontinfo;
  
  char *imagename; int ipix[2]; char label;/* read */
  int *imagetot; /* create */
  int *image;  /* alt entries */
  
  int i,j,k,fw=10,fh=10,ix[3],iy[3],ibx[3],iby[3],ifx[2],ify[2],ifxx,ifyy;
  char vstring[4][20]; int istring[3];
  double di;
  
  /* get needed stuff */
  fnames=(char **)need("p.fnames");
  fparms=(double*)need("p.fparms");
  fontinfo=(int *)find("p.fontinfo");
  if (fontinfo==0)  i=p_readfonts(0);
  fontinfo=(int *)find("p.fontinfo");
  if (arraysize("p.fontinfo")>1) {fw=fontinfo[1],fh=fontinfo[2]; }
  
  /* p.fparms: propstart propfix propend colorstart colorfix coloreend */
  /* read input */
  imagename=readname(fpin);
  ipix[0]=readint(fpin);
  ipix[1]=readint(fpin);
  label=read1charname(fpin); 
  
  printout("normal","bar image: %s size %d %d label %c for %s using p.fparms\n"
          ,imagename,ipix[0],ipix[1],label,fnames[0]);
  /* determine label strings and their sizes */
  for (i=0;i<4;i++) vstring[i][19]='\0';
  for (i=0;i<3;i+=2)
  {
    snprintf(vstring[i],18,"%lg",fparms[i]);
    snprintf(vstring[3],18,"%.3lg",fparms[i]);
    if (strlen(vstring[3])<strlen(vstring[i])) 
        snprintf(vstring[i],18,"%.3lg",fparms[i]);
    if (vstring[i][0]=='0' && vstring[i][1]!='\0') 
      for (k=1;k<strlen(vstring[i])+1;k++) vstring[i][k-1]=vstring[i][k];
    if (vstring[i][0]=='-' && vstring[i][1]=='0') 
      for (k=2;k<strlen(vstring[i])+1;k++) vstring[i][k-1]=vstring[i][k];
    istring[i]=strlen(vstring[i]);
  }
  strncpy(vstring[1],fnames[0],18);
  istring[1]=strlen(fnames[0]);
  
  /* determine location of elements and image size */
  if (label=='b' || label=='t')
  {
    j=fw*(istring[0]+istring[1]+istring[2]+2);
    if (j>ipix[0]) ipix[0]=j;
    j=fh*2.2;
    if (j>ipix[1]) ipix[1]=j;
    ix[0]=0; 
    ix[2]=ipix[0]-istring[2]*fw;
    ix[1]=.5*(fw*istring[0]+ix[2])-.5*fw*istring[1];
    ibx[0]=0; ibx[2]=ipix[0];
    if (label=='b') 
    { 
      for (i=0;i<3;i++) iy[i]=0;
      iby[0]=fh*1.2; iby[2]=ipix[1];
    }
    else /* 't' */
    {
      for (i=0;i<3;i++) iy[i]=ipix[1]-fh;
      iby[0]=0; iby[2]=ipix[1]-1.2*fh;
    }
    /* later add lrc option checks */
    /* create image */
    printout("normal","final image size: %d by %d, %s %s %s\n",ipix[0],ipix[1],vstring[0],vstring[1],vstring[2]);
    imagetot=(int*)createarray(imagename,2+ipix[0]*ipix[1],'i',0);
    image=imagetot+2;
    imagetot[0]=ipix[0]; imagetot[1]=ipix[1];
    for (i=0;i<ipix[0]*ipix[1];i++) image[i]=0;
    
    /* put on labels */
    for (i=0;i<3;i++) p_symbol(ix[i],iy[i],vstring[i],1,1,1,0,1,imagetot);
    /* put on color bar(s) */
    if (fparms[1]<=fparms[0] || fparms[1]>=fparms[2] || fparms[4]<0)
    {
      /* single bar */ 
      if (label =='b' || label=='t')
      {
        ify[0]=iby[0]; ify[1]=iby[2];        
        di=(ibx[2]-ibx[0])/(fparms[5]-fparms[3]+1);
        for (i=fparms[3];i<=fparms[5];i++) /* each color */
        { 
          ifx[0]=ibx[0]+(i-fparms[3])*di;
          if (i>fparms[3]) ifx[0]++;
          ifx[1]=ibx[0]+(i+1-fparms[3])*di;
          for (ifxx=ifx[0];ifxx<=ifx[1];ifxx++) 
            for (ifyy=ify[0];ifyy<=ify[1];ifyy++)
              p_setpixel(ifxx,ifyy,i,image,ipix);
        }
        /* later add multiple bars and lrc option */
      }
    }
    else printout("warning pc_bar","segmented color bar not yet implemented");
  }
}

