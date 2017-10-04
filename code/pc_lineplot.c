/* contains pc_lineplot   */

#include "global.h"
#include "p_plot.h"


/*-----------------------------------------------------------*/
/* this is a simple subset of an independent lineplot program
 Like the independent program, it reads a file 
 rather than using current m4d arrays - that may come later 
 note input is different ! */

void pc_lineplot(FILE *fpin, FILE *fprint)
{   
  int loop, i, j, k;   
  /* input */
  char *imagename, *datafile;
  char *comment; /* for bottom of box */ 
  double x1,x2,y1,y2; /* box info in pixels x2+20 and y2+20 are the image size */
  int icbox,icgrid; /* color box, color grid (use ticks if 0) */
  char *label[2];  double amin[2],amax[2],da[2];  /* axis info  x then y axis */
  /* repeat following input until named[0]==end  */
  char  *named[2]; double fac[2], add[2]; /* x variable then y variable  var*fac + add */
  int icolor, linew, cperline; char *symb;   /* lines and symbol */
  
  /* for data file */
  FILE *fdata;
  int nvariables, nvalues;
  char **vname, *header;
  double *values, *valuesin, *valtemp;
  int ivt1, ivt2, ivtt;
  
  /* for image*/
  int ipix[2],*imagetot,*image;
  
  /* font info */
  int fw=0, fh=0, *fontinfo;
  
  /* for box, axes, and labels */
  int ntick[2],jfont=1;  
  char vstring[20],*bottom; 
  double blog,gl;
  
  /* for graph lines */
  int iv[2];
  
  double boxll[2],boxtr[2];
  double *vplot[2] ,vmin[2],vmax[2];
  
  /*     read input names */
  imagename=readname(fpin);
  datafile=readfilename(fpin);    
  printout("normal","imagename %s datafile %s\n",imagename,datafile);
  
  /* read  datafile */
  fdata=safefopen(datafile,"r"); 
  
  header=readline(fdata);      /* comment line */
  printout("normal"," input file Header: %s\n",header);
  nvariables=readint(fdata);
  vname=(char **)tmalloca(nvariables,'p');
  printout("normal","     Variables:");
  for (k=0;k<nvariables;k++)  /*data variable names */
  {
    vname[k]=readname(fdata); 
    printout("normal"," %s",vname[k]);
  }
  ivt2=0; ivtt=nvariables*100;        /* read data in bundles of ivtt */
  valuesin=(double *)malloc(ivtt * sizeof(double));
  j=1; 
  while (j==1)
  {   
    for (ivt1=0; ivt1<ivtt; ivt1++)
    { 
      j=fscanf(fdata,"%lg",&valuesin[ivt1+ivt2]);   /* data values */
		if (j!=1) break;  
    }
    if (j==1)                     /* copy data to read more */
    {   
      ivt2=ivt1+ivt2;
      valtemp=malloc((ivt2+ivtt)*sizeof(double));
      for (ivt1=0; ivt1<ivt2; ivt1++) valtemp[ivt1]=valuesin[ivt1];
      free(valuesin); 
      valuesin=valtemp;
    }
  } /* finished reading datafile */
  fclose(fdata);
  nvalues=(ivt1+ivt2)/nvariables;
  printout("normal","    %d values\n",nvalues);
  values=malloc(nvalues*nvariables * sizeof(double));
  for  (k=0;k<nvariables;k++)  for (j=0;j<nvalues; j++)
    values[j+nvalues*k]=valuesin[k+nvariables*j];	 
  free(valuesin);
  /* finished setting datafile arrays */
  
  /*   read box info */
  comment=readname(fpin); 
  printout("normal","comment for bottom of box:  %s\n",comment);  
  x1=readdouble(fpin);
  x2=readdouble(fpin);
  y1=readdouble(fpin);
  y2=readdouble(fpin);
  icbox=readint(fpin);
  icgrid=readint(fpin);
  printout("normal","x, then y extend of box: %lg %lg %lg %lg ",x1,x2,y1,y2);
  printout("normal","boxcolor, gridcolor %d %d\n", icbox,icgrid);
  
  /* read axis info */
  for (k=0;k<2;k++)
  {  
    label[k]=readname(fpin);
    amin[k]=readdouble(fpin);
    amax[k]=readdouble(fpin);
    da[k]=readdouble(fpin);
    printout("normal","axis %d:  %s %g %g %g\n",k+1,label[k],amin[k],amax[k],da[k]);
  }
  /* get  font info */
  fontinfo=(int *)find("p.fontinfo");
  if (fontinfo==0)  i=p_readfonts(0);
  fontinfo=(int *)find("p.fontinfo");
  if (fontinfo>0) {fw=fontinfo[1],fh=fontinfo[2]; }
  
  /* create image */
  ipix[0]=x2+20;
  ipix[1]=y2+20;
  imagetot=(int*)createarray(imagename,2+ipix[0]*ipix[1],'i',0);
  image=imagetot+2;
  imagetot[0]=ipix[0]; imagetot[1]=ipix[1];
  for (i=0;i<ipix[0]*ipix[1];i++) image[i]=0;
  
  /* determine type of axes */
  for (k=0;k<2;k++)
  { 
    if (da[k]==0)  /* log scale */ 
    {  
      if (amin[k]>0)
      { blog=log10((double)amin[k]); i=blog; if (blog<0) i=blog-.999;
        amin[k]=pow(10.,(double)i);
        blog=log10((double)amax[k]);  j=blog; if (blog<0) j=blog-.999;
        amax[k]=pow(10.,(double)j);
        ntick[k]=j-i-1;
      }
      else  /* switch to regular scaling */
      {
        blog=log(amax[k]-amin[k]);  i=blog; if (blog<0) i=blog-.999;
        da[k]=pow(10.,(double)i);
      }
    }
    if (da[k]!=0)
    {if (amin[k]<0) i=amin[k]/da[k]-.999; else i=amin[k]/da[k]+.001;
      amin[k]=da[k]*i;
      if (amax[k]<0) j=amax[k]/da[k]-.001; else j=amax[k]/da[k]+.999;
      amax[k]=da[k]*j;
      ntick[k]=j-i-1;
    }
    if (k==0) printout("normal","          X-axis, %s, from %lg to %lg by %lg;",
                      label[k],amin[k],amax[k],da[k]); 
    else printout("normal","   Y-axis, %s, from %lg to %lg by %lg;\n",
                 label[k],amin[k],amax[k],da[k]);
  } 
  /* plot box and axes */
  if (icgrid>0)   /* draw grid */
  { 
    for (i=0; i<ntick[0];i++)
    { 
      gl=x1+(i+1)*(x2-x1)/((double)(ntick[0]+1));
      p_stlinex(gl,gl,y1,y2,1,icgrid,image,ipix);
    }
    for (i=0; i<ntick[1];i++)
    {  
      gl=y1+(i+1)*(y2-y1)/((double)(ntick[1]+1));
      p_stlinex(x1,x2,gl,gl,1,icgrid,image,ipix);
    }
  }
  if (icbox>0)        /* draw box */
  {
    p_stlinex(x1,x1,y1,y2,2,icbox,image,ipix);   
    p_stlinex(x1,x2,y2,y2,2,icbox,image,ipix);
    p_stlinex(x2,x2,y2,y1,2,icbox,image,ipix);  
    p_stlinex(x2,x1,y1,y1,2,icbox,image,ipix);
  }
  if (icgrid==0 && icbox>0)   /* add tick marks */
  { 
    for (i=0; i<ntick[0];i++)
    { 
      gl=x1+(i+1)*(x2-x1)/((double)(ntick[0]+1));
      p_stlinex(gl,gl,y1,y1+fw,2,icbox,image,ipix);
      p_stlinex(gl,gl,y2,y2-fw,2,icbox,image,ipix);
    }
    for (i=0; i<ntick[1];i++)
    { 
      gl=y1+(i+1)*(y2-y1)/((double)(ntick[1]+1));
      p_stlinex(x1,x1+fw,gl,gl,2,icbox,image,ipix);
      p_stlinex(x2,x2-fw,gl,gl,2,icbox,image,ipix);
    }
  }
  if (icbox>0)  /* label box  and grid  */
  {  
    for (i=0; i<ntick[0]+2;i++)
    { 
      gl=x1+i*(x2-x1)/((double)(ntick[0]+1)); /* location */
      if (da[0]==0) sprintf(vstring,"%g",	amin[0]*pow(10.,(double)i) );
      else  sprintf(vstring,"%g",amin[0]+i*da[0]);
      if (vstring[0]=='0' && vstring[1]!='\0') 
        for (k=1;k<strlen(vstring)+1;k++) vstring[k-1]=vstring[k];
      if (vstring[0]=='-' && vstring[1]=='0') 
        for (k=2;k<strlen(vstring)+1;k++) vstring[k-1]=vstring[k];
      if (i==0) p_symbolx(x1,y1,vstring,1,0,-4,0,icbox,imagetot);
      else if (i==(ntick[0]+1)) p_symbolx(x2,y1,vstring,1,-1,-4,0,icbox,imagetot);
      else  p_symbolx(gl,y1,vstring,1,0,-4,0,icbox,imagetot);		  
    }
    for (i=0; i<ntick[1]+2;i++)
    {  
      gl=y1+i*(y2-y1)/((double)(ntick[1]+1)); /* location */
      if (da[1]==0) sprintf(vstring,"%g",	amin[1]*pow(10.,(double)i) ); /* string */
      else  sprintf(vstring,"%g",amin[1]+i*da[1]);
      if (vstring[0]=='0' && vstring[1]!='\0') 
        for (k=1;k<strlen(vstring)+1;k++) vstring[k-1]=vstring[k];
      if (vstring[0]=='-' && vstring[1]=='0') 
        for (k=2;k<strlen(vstring)+1;k++) vstring[k-1]=vstring[k];
      if (i==0) p_symbolx(x1,y1,vstring,1,-2,1,0,icbox,imagetot); /* plot */
      else if (i==(ntick[1]+1)) p_symbolx(x1,y2,vstring,1,-2,-1,0,icbox,imagetot);
      else  p_symbolx(x1,gl,vstring,1,-2,0,0,icbox,imagetot);		
    }
    p_symbolx(.5*(x1+x2),y1,label[0],1,0,-6,0,icbox,imagetot); /* x-axis label */
    p_symbolx(x1,.5*(y1+y2),label[1],1,-4,0,0,icbox,imagetot); /* y-axis label */ 
    
    if (strlen(comment)*fw < x2-x1)
      p_symbolx(.5*(x1+x2),y1,comment,1,0,-8,0,icbox,imagetot);  /* put comment below box */
    else
    {
      cperline=ipix[0]/fw;
      bottom=(char *)tmalloca(cperline+2,'c');
      j=strlen(comment)/cperline; 
      if (strlen(comment)%cperline >0) j++;
      for (i=0;i<j;i++)
      {
        bottom[cperline]='\0';
        bottom[cperline+1]='\0';
        strncpy(bottom,comment+i*cperline,cperline);
        printout("normal","%d of %d lines: %s\n",i+1,j,bottom);
        p_symbolx(0,y1,bottom,1,1,-8-2*i,0,icbox,imagetot);
      }
    }
  }
  /* box and axes complete */
  
  
  /* items to be plotted  loop until variable name == end or c:
   ignore when variable not found */
  loop=100;
  while (loop)
  {
    loop--;
    for (k=0;k<2;k++) 
    {
      named[k]=readname(fpin); 
      /* check for end of plot data */
      if (strcmp(named[k],"end")==0) {loop=100; break; }
      if (strcmp(named[k],"c:")==0) { fseek(fpin,(long)(-3),1); loop=100; break; }
      fac[k]=readdouble(fpin);
      add[k]=readdouble(fpin);
    }
    if (loop==100) break;
    icolor=readint(fpin);
    linew=readint(fpin);
    symb=readname(fpin);
    printout("normal","        %s %lg %lg %s %lg %lg %d %d %s\n",named[0],add[0],fac[0],
            named[1],add[1],fac[1],icolor,linew,symb);
    
    /* find  variables to be plotted */
    iv[0]=-1; iv[1]=-1;
    for (i=0;i<nvariables;i++) 
    { 
      if (strcmp(named[0],vname[i])==0) iv[0]=i;
      if (strcmp(named[1],vname[i])==0) iv[1]=i;
    }
    if (iv[0]==-1) {printout("warning","variable %s not found\n",named[0]); continue; }
    if (iv[1]==-1) {printout("warning","variable %s not found\n",named[1]); continue; }
    
    /*data found, scale for plotting */
    vplot[0]=malloc(nvalues * sizeof(double));
    vplot[1]=malloc(nvalues * sizeof(double));
    for (k=0;k<2;k++)
    {      for (i=0;i<nvalues;i++) {vplot[k][i]=values[i+nvalues*iv[k]]*fac[k]+add[k];}
      vmin[k]=vplot[k][0];     vmax[k]=vplot[k][0];
      for (i=0;i<nvalues;i++) 
      { if (vplot[k][i]<vmin[k]) vmin[k]=vplot[k][i];
        if (vplot[k][i]>vmax[k]) vmax[k]=vplot[k][i];
      }
      printout("normal","          (%s*%g +%g) goes from %g to %g\n",vname[iv[k]],fac[k],
              add[k],vmin[k],vmax[k]);
    }
    
    /* plot data line */
    boxll[0]=x1; boxll[1]=y1; boxtr[0]=x2; boxtr[1]=y2;
    for (i=0;i<nvalues;i++)     /* for each line point replace data coor with paper coor  */
      for (k=0;k<2;k++)  
      {
        if (da[k]!=0)    /* paper coor reg plot */
          vplot[k][i]=boxll[k]+(boxtr[k]-boxll[k])*(vplot[k][i]-amin[k])/(amax[k]-amin[k]);
        else    /* log coor   if vplot<0   set to indicator, -99 */
        {
          if (vplot[k][i]>0)
            vplot[k][i]=boxll[k]+(boxtr[k]-boxll[k])*log(vplot[k][i]/amin[k])/log(amax[k]/amin[k]);
          else vplot[k][i]=-99;
        }
      }
    if (linew>0)
    {
      for (i=0;i<nvalues-1;i++) 
      { 
        if (vplot[0][i]!=-99 && vplot[1][i]!=-99 && vplot[0][i+1]!=-99 && vplot[1][i+1]!=-99)
          p_stlimited(vplot[0][i],vplot[1][i],vplot[0][i+1],vplot[1][i+1],
                      linew,icolor,x1,y1,x2,y2,imagetot);
      }
    }
    /*  plot symbols at data points */
    if (symb[0]!=' ' && symb[0]!='\0')  
      for (i=0;i<nvalues;i++)  
        if (vplot[0][i]>=x1 && vplot[0][i]<=x2 && vplot[1][i]>=y1 && vplot[1][i]<=y2)
          p_symbolx(vplot[0][i],vplot[1][i],symb,jfont,0,0,0,icolor,imagetot); 
    free(vplot[0]); free(vplot[1]);
  }
}

