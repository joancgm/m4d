/* contains p_openfile p_readfont p_symbol p_symbolx */

#include "global.h"
#include "p_plot.h"

/* for symbols on plots  
 looks for and reads font.list and listed fonts (jgm style .fmap files)
 sets up named arrays  p.fontinfo (i)   and p.fontpattern  (p = int *)
 for pixel fonts  */

/*-----------------------*/
/* open to read, file name or prefixed with names in p.defaultdir */
FILE *p_openfile(char *name)
{
  FILE *fp=NULL;
  char **defaultdir, *totname;
  int i;
  
  fp=fopen(name,"r");
  if (fp!=NULL) return fp;
  printout("normal","file %s not found, try directories in p.defaultdir\n",name);
  defaultdir=(char**)find("p.defaultdir");
  if (defaultdir==0) 
  { 
    printout("normal","string variable p.defaultdir not found, cannot find %s\n",name);
    return fp;
  }
  for (i=0;i<arraysize("p.defaultdir");i++)
  {
    if (defaultdir[i][strlen(defaultdir[i])-1]=='/')
      totname=setname2(defaultdir[i],name);
    else totname=setname3(defaultdir[i],"/",name);  
    fp=fopen(totname,"r");
    if (fp!=NULL) 
    {
      printout("normal","using file %s\n",totname);
      return fp;
    }
  }
  printout("normal","cannot find %s, in p.defaultdir directories\n",name);
  return fp;
}

/*-----------------------*/
/* read file font.list and the listed fonts */
int p_readfonts(int debug)
{   
  int numberoffonts, iwidth, iheight;
  int *fontinfo, *pat, **fontpattern;
  FILE *fp, *ff;
  char *fontname;
  int i,j,k,npat,itot;  
  char let1,let2;
  
  printout("normal","looking for file font.list\n");
  fp=p_openfile("font.list");
  if (fp==NULL) 
  {
    fontinfo=(int *)createarray("p.fontinfo",1,'i',0);
    fontinfo[0]=0;
    printout("normal","warning: there will be no letters or symbolson plots\n");
    return 0;
  }
  /* have file font.list */
  numberoffonts=readint(fp);
  if (numberoffonts<=0 || numberoffonts > 100) return 0;
  fontinfo=(int *)createarray("p.fontinfo",1+2*numberoffonts,'i',0);
  fontinfo[0]=numberoffonts;
  for (i=1;i<1+2*numberoffonts;i++) fontinfo[i]=0;
  fontpattern=(int **)createarray("p.fontpattern",127*numberoffonts+1,'p',0);

  for (i=0;i<numberoffonts;i++) 
  { 
    fontname=readname(fp);  /* read next font file name */
    printout("normal","font %d: %s\n",i+1,fontname);
    ff=p_openfile(fontname);
    if (ff==NULL) 
    {
      printout("normal"," warning, font %s omitted\n",fontname); 
      continue;
    }
    iwidth=readint(ff);
    iheight=readint(ff);
    fontinfo[1+2*i]=iwidth;
    fontinfo[2+2*i]=iheight;
    /* set up blank characters as default */
    pat=malloc(sizeof(int));
    pat[0]=iwidth*iheight;
    for (j=0;j<127;j++) fontpattern[j+127*i]=pat;
   
    /* read characters */
    for (k=0;k<127;k++)
    {
      p_newline(ff);  
      if (fscanf(ff,"%c%c",&let1,&let2)==EOF) break;
      /*  printout("normal","letter %c %c\n",let1,let2); */
      if (let1!=let2) continue;
      if ((int)let1<0 || (int)let1>=127) continue;
      npat=readint(ff);
      if (npat<=0) continue;
      pat=malloc((1+npat) * sizeof(int)); 
      itot=0;
      for (j=0;j<npat;j++) 
      {
        pat[j]=readint(ff);
        itot+=pat[j];
      }
      pat[npat]=iwidth*iheight-itot;  /* make sure there are enough items */
      fontpattern[(int)let1+127*i]=pat;
    }
  }
  printout("normal","read %d fonts\n",numberoffonts);
  return numberoffonts;
}
void p_symbolx(double x, double y,char *string,int ifont,int leftright,int updown,
               int aorp,int color, int* imagetot)  
{
  int i,j;
  i=x+.5;
  j=y+.5;
  p_symbol(i,j,string,ifont,leftright,updown,aorp,color,imagetot);
}
/* place string at pixel location x,y  */
/*.........................*/
void p_symbol(int ixin,int iyin, char *string,int ifont,int leftright,int updown,
              int aorp,int color, int* imagetot)    
{ 
  int ix,iy;
  int ist,ift, nyxset=0, nyyset=0, i,j,k, n, ipat,nzero,nycolor,*pat,fw,fh,is,ie;
  int *fontinfo, **fontpattern;
  int debug=0;
  
  ix=ixin;
  iy=iyin;
  fontinfo=(int *)find("p.fontinfo");
  fontpattern=(int**)find("p.fontpattern");
  if (fontinfo==0 || fontpattern==0) 
  {
    i=p_readfonts(debug);
    if (i==0) return;
    fontinfo=(int *)need("p.fontinfo");
    fontpattern=(int**)need("p.fontpattern");
  }
  ift=ifont-1;
  if (ift<0 || ift>fontinfo[0]-1) return;
  fw=fontinfo[1+2*ift];
  fh=fontinfo[2+2*ift];
  if (debug>0)  printout("normal","ix %d, iy %d, fw %d, fh %d\n",ix,iy,fw,fh);
  
  /* final value of ix,iy is bottom left of character */ 
  ist=strlen(string);
  /* find true center if needed */
  if (ist==1 && leftright==0 && updown==0 ) 
  {  
    int *letw, *leth;
    letw=malloc(fw*sizeof(int)); for (i=0;i<fw;i++) letw[i]=0;
    leth=malloc(fh*sizeof(int)); for (j=0;j<fh;j++) leth[j]=0;
    i=0; j=0; ipat=0; nzero=0; nycolor=0; 
    pat=fontpattern[(int)string[0]+127*ift];
    while (i<fw)
    { 
      if (pat[ipat]==0)  {nzero++;  if (nzero>2) break;}
      else
        for (n=0; n<pat[ipat];n++)
        {
          if (nycolor==1) { letw[i]=nycolor; leth[j]=nycolor;}
          j++;   if (j==fh) {j=0; i++;}
          if (i==fw) break;
        }
      ipat++;  nycolor=1-nycolor;
    }
    if (debug>1)  
    {
      printf ("center %s  ix,iy start %d %d,  footprints -\n",string,ix,iy);
      for (i=0;i<fw;i++) printout("normal","%d",letw[i]); printout("normal","\n");	  
      for (j=0;j<fh;j++) printout("normal","%d",leth[j]); printout("normal","\n");
    }
    /* leftright adjustment */
    is=fw;    for (i=0;i<fw;i++) if (letw[i]==1) {is=i; break;}
    ie=0;    for (i=fw-1; i>=0;i--) if (letw[i]==1) {ie=i; break;}
    nyxset=1;  ix=ix-fw+(ie+is+1)/2; 
    if (debug>1) printout("normal","is %d ie %d ix %d\n",is,ie,ix);
    
    /* updown adjustment */
    is=fh;  for (j=0;j<fh;j++) if (leth[j]==1) {is=j; break;}
    ie=0; 	for (j=fh-1; j>=0; j--) if (leth[j]==1) {ie=j; break;}
    nyyset=1; iy=iy-fh+(ie+is+1)/2;
    
    if (debug>1)  printout("normal","ix iy use  %d %d\n",ix,iy); 
    free (letw); free(leth);
  } /* end true character center */
  
  if (debug>0) printout("normal","ix start %d nyxset %d ",ix,nyxset);
  if (nyxset==0) /* adjust ix */
  { 
    if (leftright>1) ix=ix+(leftright-1)*fw;
    else if (leftright==0) ix=ix-.5*fw*ist;
    else if (leftright==-1) ix=ix-fw*ist;
    else  if (leftright<-1) ix=ix-fw*ist+(leftright+1)*fw;
  if (debug>0)  printout("normal","ix adjust %d,  leftright %d, fw %d, ist %d\n",ix,leftright,fw,ist);
  }
  if (debug>0) printout("normal","iy start %d nyyset %d ",iy,nyyset);
  if (nyyset==0) /* adjust iy */
  { 
    if (updown>1) iy=iy+(updown-1)*.5*fh;
    else if (updown==0) iy=iy-.5*fh;
    else if (updown==-1) iy=iy-fh;
    else if(updown<-1) iy=iy+(updown+1)*.5*fh;
  if (debug>0)	 printout("normal","iy adjust %d, updown %d, fh %d\n",iy,updown,fh);
  }
  /*  map string to image */
  for (k=0;k<ist;k++)
  { 
    pat=fontpattern[(int)string[k]+127*ift];
    i=0; j=0; ipat=0; nzero=0; nycolor=0;
    while (i<fw)
    { 
      if (pat[ipat]==0)  {nzero++;   if (nzero>2) break;  }
      else
        for (n=0; n<pat[ipat];n++)
        { 
          is=i+ix; 
          ie=iy+fh-j;
          {
            if (nycolor==1) p_setpixel(is,ie,color,imagetot+2,imagetot);
            if (nycolor==0 && aorp!=0) p_setpixel(is,ie,0,imagetot+2,imagetot);
          }
          j++;   
          if (j==fh) {j=0; i++;} 
          if (i==fw) break;
        }
      ipat++;  
      nycolor=1-nycolor;
    }
    ix=ix+fw;
  }
}


