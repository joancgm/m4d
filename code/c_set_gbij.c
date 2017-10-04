/* contains c_set_gbij */
#include "global.h"
/*  set gbij (between the points) */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_set_gbij(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep; /* needed arrays */
  double *bij;
  int *itrange; /* use if available */
  double *gbij;  /* set */
  
  int i,j,k,its,ite,ip[4],i4dp[4],i4dm[4],ii[4],mc[8],ia[3],mmid,iallm,iall;
  double cave[8],bm[3][3],bb,bbb;
  int two[3]={2,2,2};
  int ij2n[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iall=Prod4(i4d);
  
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  wherep=(int *)need("wherep");
  bij=(double *)need("bij"); 
  
  gbij=(double *)createarray("gbij",iallm,'d',0);
  
  geomcinit();
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* each cont. c.v. */
  { 
    mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3]));
    ii[3]=ip[3];        /* corner indices */
    Loop3(ia,0,two)
    { 
      Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
      mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
    }
    
    if (wherep[In4(ip,i4dp)]<0) Loop(k,0,8) cave[k]=.125;
    else geomcvave(cave,ip); /* for average values */
    
    Loop(i,0,3) Loop(j,0,3) bm[i][j]=0; 
    Loop(k,0,8)  
    Loop(i,0,3) Loop(j,0,3) bm[i][j] += cave[k]*bij[mc[k]+iall*ij2n[i][j]];
    
    bb=0; Loop(i,0,3) Loop(j,0,3) bb += bm[i][j]*bm[i][j];
    bbb=0; Loop(i,0,3) Loop(j,0,3) Loop(k,0,3) bbb += bm[i][j]*bm[j][k]*bm[k][i];
    gbij[mmid]=1.-4.5*bb+9.*bbb;
  } 	 /* each cont c.v. */
}	 
