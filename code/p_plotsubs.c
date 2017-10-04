/*   contains p_outjgm p_setpixel, p_stline p_stlinex p_stlimited p_fillx2d p_filltri2d p_painttri p_xlin */
/*   plot to a pixel image */ 
/*   the color map is specified separately, this uses color numbers */
#include "global.h"
#include "p_plot.h"
/*--------------------*/
void p_outjgm(char*filename, int*image,int*ipix)
{     
  FILE *fp;   
  int n,nbpix,ic,color,i,j;
  
  fp=safefopen(filename,"w");
  fprintf(fp,"%d %d\n",ipix[0],ipix[1]);
  
  n=0; nbpix=0; ic=image[0+ipix[0]*(ipix[1]-1)];
  for (i=0;i<ipix[0];i++) for (j=ipix[1]-1;j>=0;j--)  
  {
    color=image[i+ipix[0]*j];
    if (ic==color) nbpix++;
    else
    { fprintf(fp,"%d %d ",nbpix,ic);
      n++;   if (n>8) {fprintf(fp,"\n");  n=0;}
      ic=color; nbpix=1;
    }
    
  }
  if (nbpix>0) fprintf(fp,"%d %d",nbpix,ic);
  fprintf(fp,"\n"); 
  fclose(fp);
}
/*-----------------------*/
void p_setpixel(int i,int j,int pt_color,int *image,int *ipix)
{
  if (i>=0 && i<ipix[0] && j>=0 && j<ipix[1]) 
    image[i+ipix[0]*j]=pt_color;
}
/*------------------------*/
void p_stline(int i1,int i2,int j1,int j2,int line_width,int line_color,
              int *image,int *ipix)
{     int is,ie,js,je,i,j,k;
  /*   decide which way to go */
  is=min(i1,i2); 
  ie=i1+i2-is;  
  js=min(j1,j2); 
  je=j1+j2-js;
  if  (ie-is >= je-js) 
    for (i=is; i<=ie; i++)
    {          j=j1+(i-i1)*(j2-j1)/((float)(i2-i1));
      j-=line_width/2;   
      for (k=0;k<line_width;k++) p_setpixel(i,j+k,line_color,image,ipix);
    }
  else
    for (j=js; j<=je; j++)
    {          i=i1+(j-j1)*(i2-i1)/((float)(j2-j1));   
      i-= line_width/2;
      for (k=0;k<line_width;k++) p_setpixel(i+k,j,line_color,image,ipix);
    }
  
}
/*------------------------*/
void p_stlinex(double x1, double x2, double y1, double y2, int line_width, int line_color,
               int *image,int *ipix)
{  
  int i1,i2,j1,j2;
  int is,ie,js,je,i,j,k;
  double x,y;
  i1=x1+.5;
  i2=x2+.5;
  j1=y1+.5;
  j2=y2+.5;
  p_setpixel(i1,j1,line_color,image,ipix);
  p_setpixel(i2,j2,line_color,image,ipix);
  if (i2==i1 && j2==j1) return;
  /*   decide which way to go */
  is=min(i1,i2); 
  ie=i1+i2-is;  
  js=min(j1,j2); 
  je=j1+j2-js;
  if  (abs(x2-x1) >= abs(y2-y1)) 
    for (i=is; i<=ie; i++)
    {          
      y=y1+(i-x1)*(y2-y1)/(x2-x1)+.5;
      j=y-(line_width/2);   
      /* printout("normal"," [%d %d]",i,j); */
      for (k=0;k<line_width;k++) p_setpixel(i,j+k,line_color,image,ipix);
    }
  else
    for (j=js; j<=je; j++)
    {          
      x=x1+(j-y1)*(x2-x1)/(y2-y1)+.5;
      i=x-(line_width/2);
      /* printout("normal"," (%d %d)",i,j); */
      for (k=0;k<line_width;k++) p_setpixel(i+k,j,line_color,image,ipix);
    }
  /* printout("normal","\n"); */
}
/*-----------------------------------------------------------*/
/*calls stline for part of line in the box */
void p_stlimited(double xp1,double yp1,double xp2,double yp2,int linew, int ic, 
                double xb1, double yb1,double xb2,double yb2, int *image)  

{ 
  double fx1,fx2,fy1,fy2,x1,x2,y1,y2;
  fx1=(xp1-xb1)/(xb2-xb1);
  fx2=(xp2-xb1)/(xb2-xb1);
  fy1=(yp1-yb1)/(yb2-yb1);
  fy2=(yp2-yb1)/(yb2-yb1);
  x1=xp1; x2=xp2; y1=yp1; y2=yp2;
  if (fx1<0.) { if (fx2<0.) return;
    else           { x1=xb1;   fx1=0.; y1=yp1+(yp2-yp1)*(x1-xp1)/(xp2-xp1); } }  
  else if (fx1>1.) { if (fx2>1.) return;
    else           { x1=xb2;   fx1=1.; y1=yp1+(yp2-yp1)*(x1-xp1)/(xp2-xp1); } }    
  if (fy1<0.) { if (fy2<0.) return;
    else           { y1=yb1;   fy1=0.; x1=xp1+(xp2-xp1)*(y1-yp1)/(yp2-yp1); } } 
  else if (fy1>1.) { if (fx2>1.) return;
    else           { y1=yb2;   fy1=1.; x1=xp1+(xp2-xp1)*(y1-yp1)/(yp2-yp1); } }
  if (fx2<0.)             { x2=xb1;  y2=yp1+(yp2-yp1)*(x2-xp1)/(xp2-xp1); } 
  else if (fx2>1.)    { x2=xb2;  y2=yp1+(yp2-yp1)*(x2-xp1)/(xp2-xp1); } 
  if (fy2<0.)             { y2=yb1;  x2=xp1+(xp2-xp1)*(y2-yp1)/(yp2-yp1); } 
  else if (fy2>1.)     { y2=yb2;  x2=xp1+(xp2-xp1)*(y2-yp1)/(yp2-yp1); } 
  p_stlinex(x1,x2,y1,y2,linew,ic,image+2,image);
  return;
}

/*-----------------------------------------------------------*/
/* fill the triangle with color ic, note y1=y2=y12 */
void p_fillx2d(double x1, double x2, double x3, double y12, double y3,
               int ic, int *image,int *ipix)
{  
  int i,j,js,je,is,ie;
  double ymin,ymax,xmin,xmax,xa,xb,x,y;
  
  if (y12 > y3) { ymin=y3; ymax=y12; }
  else          { ymin=y12; ymax=y3; }
  js=ymin; if ( (double)js<ymin ) js=js+1;
  je=ymax; if ( (double)je>ymax ) je=je-1; 
  for (j=js; j<je+1; j++)
  { 
    y=j;
    xa=p_xlin(x1,y12,x3,y3,y);
    xb=p_xlin(x2,y12,x3,y3,y);
    if (xb > xa) { xmin=xa; xmax=xb; }
    else         { xmin=xb; xmax=xa; }
    is=xmin; if ( (double)is<xmin ) is=is+1;
    ie=xmax; if ( (double)ie>xmax ) ie=ie-1;
    for(i=is; i<ie+1; i++)
    {  x=i;
      p_setpixel(i,j,ic,image,ipix); 
    } 
  }
}
/*-----------------------------------------------------------*/
/* break triangle up into sub triangles, two points with same y */ 
void p_filltri2d(double x1, double x2, double x3, double y1, double y2, double y3,
                 int ic, int *image,int *ipix)
{   
  double ymax,ymin,x; int imin,imax;
  
  ymin=y1; imin=1; 
  if (y2<ymin) {ymin=y2; imin=2; }
  if (y3<ymin) {ymin=y3; imin=3; }
  ymax=y3; imax=3;
  if (y2>ymax) {ymax=y2; imax=2; }
  if (y1>ymax) {ymax=y1; imax=1; }
  if (imin!=1 && imax!=1)
  { x=p_xlin(x2,y2,x3,y3,y1);
    p_fillx2d(x1,x,x3,y1,y3,ic,image,ipix);
    p_fillx2d(x1,x,x2,y1,y2,ic,image,ipix);
  }
  else if (imin!=2 && imax!=2)
  { x=p_xlin(x1,y1,x3,y3,y2);
    p_fillx2d(x2,x,x3,y2,y3,ic,image,ipix);
    p_fillx2d(x2,x,x1,y2,y1,ic,image,ipix);
  }
  else 
  { x=p_xlin(x1,y1,x2,y2,y3);
    p_fillx2d(x3,x,x2,y3,y2,ic,image,ipix);
    p_fillx2d(x3,x,x1,y3,y1,ic,image,ipix);
  }
}
/*-----------------------------------------------------------*/
/* linearly interpolate for x, given y and two points on the line */
double p_xlin(double x1, double y1, double x2, double y2, double y)
{ 
  double x;
  if (y1==y2) return x1;
  else
  { x=x1+(y-y1)*(x2-x1)/(y2-y1);
    return x;
  }
}
/*-----------------------------------------------------------*/
/* break a paint triangle up into individual colors segments */
void p_painttri(double *xys[3],double pp[3],char nyout, 
              int* ic, double *pc, int ncolor, int *image, int ipix[2])
{
  int i,k,m[3],is,knum;
  double xy[3][2],p[3],xypl[5][2],f,f2;
  
  /* order points */
  m[0]=0; m[2]=2;
  for (i=0;i<3;i++) 
  {
    if (pp[i]<pp[m[0]]) m[0]=i;
    if (pp[i]>pp[m[2]]) m[2]=i;
  }
  m[1]=3-m[0]-m[2];
  
  
  for (i=0;i<3;i++)
  {
    xy[i][0]=xys[m[i]][0];
    xy[i][1]=xys[m[i]][1];
    p[i]=pp[m[i]];
  }
  /* printout("normal"," xyz %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",xy[0][0],xy[1][0],xy[2][0],
         xy[0][1],xy[1][1],xy[2][1],p[0],p[1],p[2]); */
  
  
  for (k=0;k<ncolor;k++) /* for each color fill */
  {
    if (ic[k]==0) continue; /* leave 0 'transparent' */
    
    /* initialize contour limits */
    for (i=0;i<2;i++)
    {
      xypl[0][i]=xy[0][i];
      xypl[1][i]=xy[1][i];
      xypl[2][i]=xy[2][i];
      xypl[3][i]=xy[2][i];
      xypl[4][i]=xy[0][i];
    }
    /* reset lower limit */
    
    if (k>0 || nyout=='n')
    {
      if (p[2]<= pc[k]) continue;
      if (p[0] < pc[k])
      {
        f=(pc[k]-p[0])/(p[2]-p[0]);
        for (i=0;i<2;i++) xypl[4][i]=xy[0][i]+f*(xy[2][i]-xy[0][i]);
        is=0;                                               
        if (p[1]<pc[k]) is=1;
        f2=(pc[k]-p[is])/(p[is+1]-p[is]);
        for (i=0;i<2;i++) 
        {
          xypl[0][i]=xy[is][i]+f2*(xy[is+1][i]-xy[is][i]);
          xypl[is][i]=xypl[0][i];
        }
        /* printout("normal","low k %d f is f2 %lg %d %lg\n",k,f,is,f2); */
      }
    }
    
    /* reset upper limit */
    if (k<ncolor-1 || nyout=='n')
    {
      if (p[0]>=pc[k+1]) continue;
      if (p[2] > pc[k+1])
      {
        f=(pc[k+1]-p[0])/(p[2]-p[0]); 
        for (i=0;i<2;i++) xypl[3][i]=xy[0][i]+f*(xy[2][i]-xy[0][i]);
        is=2;  
        if (p[1]>pc[k+1]) is=1;
        f2=(pc[k+1]-p[is-1])/(p[is]-p[is-1]);
        for (i=0;i<2;i++) 
        {
          xypl[2][i]=xy[is-1][i]+f2*(xy[is][i]-xy[is-1][i]);
          xypl[is][i]=xypl[2][i];
        }
         /* printout("normal","hi k %d f is f2 %lg %d %lg\n",k,f,is,f2); */
      }
    }
    /* for (i=0;i<5;i++) printout("normal","  %lg %lg",xypl[i][0],xypl[i][1]); */
    /* eliminate duplicate points */
    knum=0;
    for (i=1;i<5;i++)
    {
      if (xypl[i][0]==xypl[knum][0] && xypl[i][1]==xypl[knum][1]) continue;
      knum++;
      xypl[knum][0]=xypl[i][0]; xypl[knum][1]=xypl[i][1]; 
    }
    if (xypl[0][0]==xypl[knum][0] && xypl[0][1]==xypl[knum][1]) knum--;
   /*  printout("normal","knum %d %d\n",knum,k); */
    if (knum<2) continue;
    for (i=2;i<=knum;i++) p_filltri2d(xypl[0][0],xypl[i-1][0],xypl[i][0],
                                      xypl[0][1], xypl[i-1][1], xypl[i][1],ic[k],image,ipix);
  }
}
