/* contains c_aveijk */
#include "global.h"
/* aveijk  namefrom nameto (on or between the points arrays)  sym[3] (i,j,k)
 n,s,m,M  n-none, s-sum pt-average, a-abc average, m-mirror + ave, M-mirror - ave
 e.g. c: aveijk U1 U1smm s m m U2 U2sMm s M m U3 U3smM s M m ""    for square duct 
 warning for mirror options it will be assumed that the grid spacing mirrors !!!
 */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_aveijk(FILE *fpin, FILE *fprint)
{ 
  int *i4d; double *a[3]; /* needed arrays */
  char *namefrom, *nameto, sym[3]; /* read */
  
  int i4dm[4],iall,iallm,*i4dx;
  int i,j,n,ii[4],jj[4],kk[4],iptf,iptt,ia[3],iamax[3],size;
  double *vf, *vto, da, *atot, vpt;
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iallm=Prod4(i4dm);
  iall=Prod4(i4d);
  
  a[0]=(double *)need("abcd");
  a[1]=a[0]+i4d[0];
  a[2]=a[1]+i4d[1];
  atot=(double *)tmalloca(iall,'d');
  
  n=20;
  while(n)
  {  
    namefrom=readname(fpin);
    if (namefrom[0]=='\0')  break;
    if (strcmp(namefrom,"c:")==0) break;
    nameto=readname(fpin);
    Loop(i,0,3) sym[i]=read1charname(fpin);
    vf=(double *)find(namefrom);
    if (vf==0) 
    {
      printout("warning c_aveijk"," %s not found, no averaging\n",namefrom);
      n--; continue;
    }
    size=arraysize(namefrom);
    if (size != iall && size !=iallm)
    {	
      printout("warning c_aveijk"," %s not on or between pts array, no averaging\n",namefrom);
      continue;
    }
    printout("normal","average %s to %s using %c %c %c\n",
            namefrom,nameto,sym[0],sym[1],sym[2]);
    vto=(double*)createarray(nameto,size,'d',0);
    Loop(i,0,3) {iamax[i]=1; if (sym[i]=='m' || sym[i]=='M') iamax[i]=2;}
    if (size==iall) i4dx=i4d; else i4dx=i4dm;
    
    /* do sums */
    Loop(i,0,size) {vto[i]=0; atot[i]=0;}
    
    Loop(ii[3],0,i4dx[3]) Loop3(ii,0,i4dx)
    Loop3(ia,0,iamax)
    {  
      jj[3]=ii[3]; kk[3]=0;
      Loop(i,0,3) {jj[i]=ii[i]; if (ia[i]==1) jj[i]=i4dx[i]-1-ii[i];}
      Loop(i,0,3) if (sym[i]=='I' && i4dx[i]==i4dx[(i+1)%3])
      {j=jj[i]; jj[i]=jj[(i+1)%3]; jj[(i+1)%3]=j; }
      iptf=In4(jj,i4dx); vpt=vf[iptf];
      Loop(i,0,3) if (ia[i]==1 && sym[i]=='M') vpt=-vpt;
      Loop(i,0,3) {kk[i]=ii[i]; if (sym[i]=='s' || sym[i]=='a') kk[i]=0; }
      iptt=In4(kk,i4dx);
      da=1;
      Loop(i,0,3) if (sym[i]=='a')
      { 
        if (size==iall) da*= a[i][min(ii[i]+1,i4d[i]-1)]-a[i][max(ii[i]-1,0)];
        else da*= a[i][ii[i]+1]-a[i][ii[i]];
      }
      atot[iptt]+=da;
      vto[iptt]+=da*vpt;
    }
    /* averages */
    Loop(ii[3],0,i4dx[3]) Loop3(ii,0,i4dx)
    {  
      iptt=In4(ii,i4dx);
      if (atot[iptt]>0) vto[iptt]/=atot[iptt];
      else   /* summed direction copy from 0 */
      { 
        kk[3]=0;
        Loop(i,0,3) {kk[i]=ii[i]; if (sym[i]=='s' || sym[i]=='a') kk[i]=0; }
        vto[iptt]=vto[In4(kk,i4dx)];
      }
    }
  }
  if (n==0) printout("warning c_aveijk"," 20 names not found, aveijk abandoned\n");
}