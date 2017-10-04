/* contains pc_gridinfo */
#include "global.h"
/*  create idim4ddouble abcddouble cltdouble to use for plotting double grid  */

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])

/* ---------------------------- */
void pc_gridinfo(FILE *fpin, FILE *fprint)
{
  int *i4d; double *a[4]; char *clt; /* need */
  int *idim4dd; double *adouble[4]; char *cltdouble; /* set */
  
  int i,j,ii[4],ialld,iid[4],iter,ich,iunset,icen,iplus,iminus;
   
  /* get needed arrays */
  i4d=(int *)need("idim4d");
  a[0]=(double *)need("abcd");
  for (i=1;i<4;i++) a[i]=a[i-1]+i4d[i-1];
  clt=(char *)need("clt");
  
  idim4dd=(int *)createarray("idim4ddouble",4,'i',0);
  for (i=0;i<3;i++) idim4dd[i]=2*i4d[i]-1;
  idim4dd[3]=1;
  
  i=idim4dd[0]+idim4dd[1]+idim4dd[2]+idim4dd[3];
  adouble[0]=(double *)createarray("abcddouble",i,'d',0);
  for (i=1;i<4;i++) adouble[i]=adouble[i-1]+idim4dd[i-1];
  adouble[3][0]=a[3][0];
  for (i=0;i<3;i++)
  {
    for (j=0;j<i4d[i];j++) adouble[i][2*j]=a[i][j];
    for (j=1;j<i4d[i];j++) adouble[i][2*j-1]=.5*(a[i][j-1]+a[i][j]);
  }
  
  ialld=Prod4(idim4dd);
  cltdouble=(char *) createarray("cltdouble",ialld,'c',0);
  for (i=0;i<ialld;i++) cltdouble[i]='*';
  ii[3]=0; iid[3]=0;
  /* same points */
 Loop3(ii,0,i4d) 
  {
    for (i=0;i<3;i++) iid[i]=2*ii[i];
    cltdouble[In4(iid,idim4dd)]=clt[In4(ii,i4d)];
  } 
  /* between points */
  for (iter=0; iter<5; iter++)
  {
    ich=0; iunset=0;
    Loop3(iid,0,idim4dd)
    {
      icen=In4(iid,idim4dd);
      if (cltdouble[icen] != '*') continue;
      iunset++;
      for (i=0;i<3;i++)
      {
        if (iid[i]==0 || iid[i]==idim4dd[i]-1) continue;
        iid[i]--; iminus=In4(iid,idim4dd); iid[i]++;
        iid[i]++; iplus=In4(iid,idim4dd); iid[i]--;
        if (cltdouble[iminus]=='*' || cltdouble[iplus]=='*') continue;
        
        if (cltdouble[iminus]==cltdouble[iplus]) 
          cltdouble[icen]=cltdouble[iminus];
        else if (cltdouble[iminus]=='s' || cltdouble[iplus]=='s') 
          cltdouble[icen]='s';
        else  cltdouble[icen]='f'; 
        ich++;
        iunset--; 
        break; 
      }
    }
    if (ich==0 || iunset==0) break;
  }
    printout("warning from pc_gridinfo","some cltdouble values not setich %d, iunset %d\n",ich,iunset);
  /* reset grid points to x */
  Loop3(ii,0,i4d) 
  {
    for (i=0;i<3;i++) iid[i]=2*ii[i];
    cltdouble[In4(iid,idim4dd)]='x';
  } 
}

  