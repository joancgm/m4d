/* contains valueatp c_valueat */
#include "global.h"
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])

/* find value in a type double array -  on, mid, or p-points */
/* nameout = namearray at specified abc  using lin interp in a, b and c */

/*--------------- return value at p   given ii,ff for on points */
double valueatp(double *v,double abc[3],double **a,int i4d[4],
                double *cvdc,int i4dp[4],int ii[3],double ff[3])
{
  int i,iall;
  double fc,flow,fhi,amid,alow,ahi;
  
  iall=Prod4(i4d);
  
  for (i=0;i<3;i++)
  {
    /* first check out of range */
    if (abc[i]<=a[i][0]) continue;
    else if (abc[i]>=a[i][i4d[i]-1]) ii[i]++;
    else  /* need to interpolate */
    {
      fc=.5; flow=.5; fhi=.5;
      if (cvdc>0) /* reset fc,flow and fhi */
      {
        fc=cvdc[In4(ii,i4d)+iall*i];
        if (ii[i]>0) {ii[i]--; flow=cvdc[In4(ii,i4d)+iall*i]; ii[i]++; }
        else if (ii[i]<i4d[i]-1) {ii[i]++; fhi=cvdc[In4(ii,i4d)+iall*i]; ii[i]--; }
      }  
      /* determine a values for fc flow and fhi */
      amid=a[i][ii[i]]+fc*(a[i][ii[i]+1]-a[i][ii[i]]);
      alow=a[i][0]; 
      if (ii[i]>0) alow=a[i][ii[i]-1]+flow*(a[i][ii[i]]-a[i][ii[i]-1]);
      ahi=a[i][i4d[i]-1];
      if (ii[i]<i4d[i]-2) ahi=a[i][ii[i]+1]+fhi*(a[i][ii[i]+2]-a[i][ii[i]+1]);
      /* now assuming space and abc interp the same reset ff and ii for p array */
      if (abc[i]<=amid) 
        ff[i]=(abc[i]-alow)/(amid-alow);
      else
      {
        ii[i]++;
        ff[i]=(abc[i]-amid)/(ahi-amid);
      }
    }
  }
  return(Value3d(v,i4dp,ii,ff)); 
}
/* -------------------------------*/
void c_valueat(FILE *fpin, FILE *fprint)
{   
  char *nameout,*namearray;  /* input */
  double abc[3];
  /* nameout[0]=namearray[ at abc] */
  
  int *i4d; double *a[3], *v; /* needed */
  double *cvdc=0; /* use for p-point interp if available */
  double *vout;  /* created */
  
  int size,iall,i4dm[4],iallm,i4dp[4],iallp,i,ii[4];
  double ff[3];
  
  nameout=readname(fpin);
  namearray=readname(fpin);
  for (i=0;i<3;i++) abc[i]=readdouble(fpin);
  printout("normal","%s[1] = %s[at %lg,%lg,%lg]\n",
          nameout,namearray,abc[0],abc[1],abc[2]);
  v=(double *)need(namearray);
  size=arraysize(namearray);
  vout=(double *)createarray(nameout,1,'d',0);
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double *)need("abcd"); 
  a[1]=a[0]+i4d[0];
  a[2]=a[1]+i4d[1];
  
  for (i=0;i<3;i++) ii[i]=findex(abc[i],a[i],i4d[i],&ff[i]);
  /*  ---------on points */
  if (size==iall)  
  {
    vout[0]=Value3d(v,i4d,ii,ff);
    return;
  }
  /* ----------- mid points */
  for (i=0;i<3;i++) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  if (size==iallm) 
  {
    for (i=0;i<3;i++)
    {
      if (abc[i]<=.5*(a[i][0]+a[i][1])) {ii[i]=0; ff[i]=0;}
      else if ( abc[i]>=.5* (a[i][i4d[i]-2] + a[i][i4d[i]-1] ) )
      {ii[i]=i4d[i]-3; ff[i]=1;}
      else 
      { 
        if (abc[i]<.5*(a[i][ii[i]]+a[i][ii[i]+1])) ii[i]--;
        ff[i]=(abc[i]-.5*(a[i][ii[i]]+a[i][ii[i]+1]))
        /(.5*(a[i][ii[i]+2]-a[i][ii[i]]));
      }    
    }
    vout[0]=Value3d(v,i4dm,ii,ff); 
    return;
  }
  /* ------------- p-points */
  for (i=0;i<3;i++) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  if (size==iallp)  
  {
    cvdc=(double *)find("cvdc");
    vout[0]=valueatp(v,abc,a,i4d,cvdc,i4dp,ii,ff);
    return;
  }
  /* not on, mid or p points */
  printout("error c_valueat","error - array %s not on,mid or p-points array\n",namearray);
  exitm4d(0);
}
