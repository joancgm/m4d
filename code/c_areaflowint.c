/* contains areaflowint, c_areaflowint */
#include "global.h"
/* area integral, L direction over irange , note simple 4-pt averages are used  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) {*(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b);} 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

int areaflowint(int L, int (*irange)[2], int *i4d, double **x, double *rho, 
                double **u, double **prop, int iprop, double *aves)
{ int L2, L3,i,j,ii[4],ix[4],ia,ib,ipt,ip;
  double ru[3],vav[20],dx1[3],dx2[3],area[3],add,flow;
  
  L2=(L+1)%3;
  L3=(L+2)%3;
  ii[3]=irange[3][0]; ix[3]=ii[3]; 
  ii[L]=irange[L][0]; ix[L]=ii[L];
  Loop(i,0,7*(iprop+1)) aves[i]=0;
  
  Loop(ii[L2],irange[L2][0],irange[L2][1])
  Loop(ii[L3],irange[L3][0],irange[L3][1])
  { 
    Loop(i,0,iprop) vav[i]=0;
    Loop(i,0,3) {dx1[i]=0; dx2[i]=0; ru[i]=0; }
    Loop(ia,0,2) Loop(ib,0,2)
    { 
      ix[L2]=ii[L2]+ia;
      ix[L3]=ii[L3]+ib;
      ipt=In4(ix,i4d);
      Loop(i,0,iprop) vav[i]+=.25*prop[i][ipt];
      Loop(i,0,3) 
      { 
        dx1[i]+=.5*(2*ia-1)*x[i][ipt];
        dx2[i]+=.5*(2*ib-1)*x[i][ipt];
        ru[i]+=.25*rho[ipt]*u[i][ipt];
      }
    } /* ia,ib */	     
    Cross(dx1,dx2,area);
    
    add=sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
    aves[0]+=add;
    Loop(j,0,3) aves[j+1]+=area[j];
    Loop(i,0,iprop) 
    { 
      aves[7*(i+1)]+=add*vav[i];
      Loop(j,0,3) aves[j+1+7*(i+1)]+=area[j]*vav[i];
    }
    flow=Dot(area,ru);
    ip=4; if (flow<0) ip=5;
    aves[ip]+=flow;
    Loop(i,0,iprop) aves[ip+7*(i+1)]+=flow*vav[i];
  }
  Loop(i,0,iprop) aves[6+7*(i+1)]=aves[4+7*(i+1)]+aves[5+7*(i+1)];
  aves[6]=aves[4]+aves[5];
  Loop(j,0,7) if (aves[j]!=0) 
    Loop(i,0,iprop) aves[j+7*(i+1)]/=aves[j];
  if (aves[0]==0) return 0; else return 1;
} 
/* perform a flow and area averages of on-the-grid-points properties
 over the grid plane and range specified */
/* input abcijk (character, plane type) 
 if a b or c  abcvalue, bcastart, bcaend, cabstart, cabend
 else if i j or k  ijkvalue, jkistart, jkiend, kijstart, kjiend 
 propname
 repeat propname until propname="" or propname not found
 max read 20 propnames
 print abcvalue or ijkvalue, area flow+ flow- flow ITER TIME (if available)
 then for each propname 
 abcvalue or ijkvalue, areaint flow+int flow-int flowint ITER TIME (if available)
 Note a b and c values will be corrected to the nearest grid values if no identical match
 keywords AFINT and the property names will also be printed (for catching with grep)
 also creates the array averages[7,21] 
 7 values overall and then for each prop are: area areax areay areaz flow+ flow- flow
 */
void c_areaflowint(FILE *fpin, FILE *fprint)
{ int *i4d; double *a[4],*x[3],*rho,*u[3]; /* need */
  int *ITER, *itrange; double *TIME; /* use if available */
  double *averages; /* create */
  
  char abcijk;  double arange[3][2];  int irange[4][2]; /* input parameters */
  char *name[20];
  
  char aori,abc[3]={'a','b','c'},ijk[3]={'i','j','k'};
  int L,L2,L3,iall,iname,i,j;
  double *prop[20],f[1];
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double *)need("abcd");
  a[1]=a[0]+i4d[0]; a[2]=a[1]+i4d[1];
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[1]+iall;
  rho=(double *)need("rho");
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  itrange=(int *)find("itrange");
  irange[3][0]=i4d[3]-1;
  if (itrange>0) irange[3][0]=itrange[1];
  ITER=(int *)find("ITER");
  TIME=(double *)find("TIME");
  averages=(double *)createarray("averages",7*21,'d',0);
  Loop(i,0,7*21) averages[i]=0;
  
  abcijk=read1charname(fpin);    /* read input */
  if (abcijk=='a') {L=0; aori='a';}
  else if (abcijk=='b') {L=1; aori='a';}
  else if (abcijk=='c') {L=2; aori='a';}
  else if (abcijk=='i') {L=0; aori='i';}
  else if (abcijk=='j') {L=1; aori='i';}
  else if (abcijk=='k') {L=2; aori='i';}
  else {printout("warning c_areaflowing","input error abcijk= %c, return\n",abcijk); return; }
  L2=(L+1)%3;
  L3=(L+2)%3;
  if (aori=='a') 
  {
    arange[L][0]=readdouble(fpin); 
	 arange[L][1]=arange[L][0];
	 arange[L2][0]=readdouble(fpin);
	 arange[L2][1]=readdouble(fpin);
	 arange[L3][0]=readdouble(fpin);
	 arange[L3][1]=readdouble(fpin);
	 Loop(i,0,3) Loop(j,0,2) irange[i][j]=findex(arange[i][j],a[i],i4d[i],f)+f[0]+.5;
	 printout("normal","integration for %c=%lg, %c=%lg to %lg, %c=%lg to %lg\n",
            abc[L],a[L][irange[L][0]],abc[L2],a[L2][irange[L2][0]],a[L2][irange[L2][1]],
            abc[L3],a[L3][irange[L3][0]],a[L3][irange[L3][1]]);
  }
  else
  {
    irange[L][0]=readint(fpin);
	 irange[L][1]=irange[L][0];
	 irange[L2][0]=readint(fpin);
	 irange[L2][1]=readint(fpin);
	 irange[L3][0]=readint(fpin);
	 irange[L3][1]=readint(fpin);
	 Loop(i,0,3) Loop(j,0,2) irange[i][j]=max(0,min(i4d[i]-1,irange[i][j]));
	 printout("normal","integration for %c=%d, %c=%d to %d, %c=%d to %d\n",
            ijk[L],irange[i][0],ijk[L2],irange[L2][0],irange[L2][0],
            ijk[L3],irange[L3][0],irange[L3][0]);
  }
  printout("normal","properties:");
  Loop(iname,0,20)   /* read and find properties */
  { 
    name[iname]=readname(fpin);
    prop[iname]=(double *)find(name[iname]);
    if (prop[iname]==0)  break;
    printout("normal"," %s",name[iname]);
  }
  printout("normal","\n");
  
  areaflowint(L,irange,i4d,x,rho,u,prop,iname,averages);
  
  printout("normal","   a,b,or c  area  flow+ flow- flow  (property averages)\n");
  Loop(j,0,iname+1)  /* print results */
  { 
    printout("normal","AFINT");
    if (j==0) printout("normal"," %c%lg",abc[L],a[L][irange[L][0]]);
    else  printout("normal"," %s_%c%lg",name[j-1],abc[L],a[L][irange[L][0]]);
    printout("normal"," %lg  %lg %lg  %lg",
            averages[7*j],averages[4+7*j],averages[5+7*j],averages[6+7*j]);
    if (ITER>0) printout("normal"," ITER %d",ITER[0]);
    if (TIME>0) printout("normal"," TIME %lg",TIME[0]);
    printout("normal","\n");
  }
}