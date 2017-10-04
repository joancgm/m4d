/*   contains p_fillcontour */
#include "global.h"
#include "p_plot.h"
/* fill type contours */
/* p.fnames[3]: prop name, clt types to be filled (i-pt pass), 
 nyout - n/y fill out of range with end color
 p.fparms[6]: props propf prope colors colorf colore
 propf is at color f unless colorf is negative, or propf is out of range
 contours from props to prope, prope >= props
 if they are both equal, they are set to cover the range
 and the values in p.fparms are reset. */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Vabs(a) (sqrt((double)( *(a) * *(a) + *((a)+1) * *((a)+1) + *((a)+2) * *((a)+2) )))

void p_fillcontour(double *p, int nystep, int iadd, char ** fnames, double *fparms,
                   int*i4d, double *x[3],char *clt,
                   int irange[3][2],double *gxx[2],double scale,double offset[2],
                   int *image, int ipix[2])
{
  int ii[4]={0,0,0,0},jj[4]={0,0,0,0},L,L2,L3,i,j,k,i4dc[4]={1,1,1,1},mm[4];
  double ps,pf,pe,colorfd,dp1,dp2;
  int colors,colorf,colore,ncolor,ncolor1,ncolor2,ncltdo,nydo;
  char nyout,*cltdo;
  int *ic; double *pc;
  double xpl[5][2],xpt[4][3],pvalue[5],*xy[3],ptri[3];
  
  if (nystep==0) for (i=0;i<3;i++) i4dc[i]=i4d[i];
  else for (i=0;i<3;i++) i4dc[i]=i4d[i]-1+2*iadd;
  
  printout("normal","i4dc %d %d %d %d\n",i4dc[0],i4dc[1],i4dc[2],i4dc[3]);
  ps=fparms[0];
  pf=fparms[1];
  pe=fparms[2];
  colors=fparms[3]+.01;
  colorf=fparms[4]+.01;
  colorfd=fparms[4];
  colore=fparms[5]+.01;
  cltdo=fnames[1];
  ncltdo=strlen(cltdo);
  nyout=fnames[2][0];
  /* add clt types */
  
  /* choose direction */
  L=0;
  for (i=1;i<3;i++) if (irange[i][1]-irange[i][0] < irange[L][1]-irange[L][0]) L=i;
  L2=(L+1)%3;
  L3=(L+2)%3;
  k=0;
  if (ps==pe) /* autoscale */
  {
    for (ii[L]=irange[L][0];ii[L]<=irange[L][1];ii[L]++)
      for (ii[L2]=irange[L2][0];ii[L2]<irange[L2][1];ii[L2]++)
        for (ii[L3]=irange[L3][0];ii[L3]<irange[L3][1];ii[L3]++) 
        {
          /*  printout("normal","ijk %d %d %d\n",ii[0],ii[1],ii[2]); */
          for (i=0;i<3;i++) jj[i]=ii[i];
          mm[0]=In4(jj,i4d);
          jj[L2]++; mm[1]=In4(jj,i4d);
          jj[L3]++; mm[2]=In4(jj,i4d);
          jj[L2]--; mm[3]=In4(jj,i4d);
          jj[L3]--;
          nydo=0;
          for (j=0;j<ncltdo;j++) for (i=0;i<4;i++) if (clt[mm[i]]==cltdo[j]) nydo=1;
          if (nydo==0) continue;
          
          if (nystep==0) for (j=0;j<4;j++) 
          {
            if (k==0) { ps=p[mm[j]]; pe=p[mm[j]]; k=1;}
            else {ps=min(ps,p[mm[j]]); pe=max(pe,p[mm[j]]); }
          }
          else
          {
            for (i=0;i<3;i++) jj[i]+=iadd;
            pvalue[4]=p[In4(jj,i4dc)];
            if (k==0) { ps=pvalue[4]; pe=pvalue[4]; k=1;}
            else {ps=min(ps,pvalue[4]); pe=max(pe,pvalue[4]); }
          }
        }
    printout("normal"," autoscale contours from %lg to %lg\n",ps,pe);
    if (ps==pe)
    {
      printout("normal","no range in the property, contours omitted\n");
      return;
    }
    fparms[0]=ps;
    fparms[2]=pe;
  }
  /* determine color and contour list */
  if (pf<=ps || pf>=pe || colorfd<0) 
  {
    ncolor=abs(colore-colors)+1;
    dp1=ncolor; dp1=(pe-ps)/dp1;
    ic=tmalloca(ncolor,'i');
    k=1; if (colore<colors) k=-1;
    for (i=0;i<ncolor;i++) ic[i]=colors+i*k;
    pc=tmalloca(ncolor+1,'d');
    for (i=0;i<ncolor+1;i++) pc[i]=ps+dp1*i;
  }
  else /* segmented color range */
  {
    ncolor1=abs(colors-colorf)+1;
    ncolor2=abs(colore-colorf)+1;
    ncolor=ncolor1+ncolor2;
    ic=tmalloca(ncolor,'i');
    k=1; if (colorf<colors) k=-1;
    for (i=0;i<ncolor1;i++) ic[i]=colors+i*k;
    k=1; if (colore<colorf) k=-1;
    for (i=0;i<ncolor2;i++) ic[i+ncolor1]=colorf+i*k;
    pc=tmalloca(ncolor+1,'d');
    dp1=ncolor1; dp1=(pf-ps)/dp1;
    for (i=0;i<ncolor1+1;i++) pc[i]=ps+dp1*i;
    dp2=ncolor2; dp2=(pe-pf)/dp2;
    for (i=0;i<ncolor2+1;i++) pc[i+ncolor1]=pf+dp2*i;
  }
  for (i=0;i<ncolor;i++) printout("normal","%lg %d  ",pc[i],ic[i]);
   printout("normal","%lg\n",pc[ncolor]); 
  
  /* ready to do contours */
  for (ii[L]=irange[L][0];ii[L]<=irange[L][1];ii[L]++)
    for (ii[L2]=irange[L2][0];ii[L2]<irange[L2][1];ii[L2]++)
      for (ii[L3]=irange[L3][0];ii[L3]<irange[L3][1];ii[L3]++) 
      {
        /*  printout("normal","ijk %d %d %d\n",ii[0],ii[1],ii[2]); */
        for (i=0;i<3;i++) jj[i]=ii[i];
        mm[0]=In4(jj,i4d);
        jj[L2]++; mm[1]=In4(jj,i4d);
        jj[L3]++; mm[2]=In4(jj,i4d);
        jj[L2]--; mm[3]=In4(jj,i4d);
        jj[L3]--;
        nydo=0;
        for (j=0;j<ncltdo;j++) for (i=0;i<4;i++) if (clt[mm[i]]==cltdo[j]) nydo=1;
        if (nydo==0) continue;
        xpl[4][0]=0; xpl[4][1]=0;
        for (j=0;j<4;j++)
        {
          for (i=0;i<3;i++) xpt[j][i]=x[i][mm[j]];
          for (i=0;i<2;i++) xpl[j][i]=(Dot(xpt[j],gxx[i])+offset[i])*scale;
          for (i=0;i<2;i++) xpl[4][i]+=.25*xpl[j][i];
        }
        pvalue[4]=0;
        if (nystep==0) for (j=0;j<4;j++) { pvalue[j]=p[mm[j]]; pvalue[4]+=.25*pvalue[j]; }
        else 
        { 
          for (i=0;i<3;i++) jj[i]+=iadd;
          pvalue[4]=p[In4(jj,i4dc)];
          for (j=0;j<4;j++) pvalue[j]=pvalue[4];         
        }
        for (i=0;i<4;i++) 
        {
          xy[0]=xpl[i]; xy[1]=xpl[(i+1)%4]; xy[2]=xpl[4]; 
          ptri[0]=pvalue[i]; ptri[1]=pvalue[(i+1)%4]; ptri[2]=pvalue[4];
          /*  printout("normal"," call paintr %lg %lg %lg \n",ptri[0],ptri[1],ptri[2]); */
          p_painttri(xy,ptri,nyout,ic,pc,ncolor,image,ipix);
        }
      }
}


