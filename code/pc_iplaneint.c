/* contains pc_iplaneint */
#include "global.h"
/* area integral, L direction over irange , note simple 4-pt averages are used  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) {*(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b);} 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

int areaflowint(int L, int (*irange)[2], int *i4d, double **x, double *rho, 
                double **u, double **prop, int iprop, double *aves);

/* array averages[7,21] 
 7 values overall and then for each prop are: area areax areay areaz flow+ flow- flow
 */
/* area and mass average properties over i-planes and send to lineplot file */
/* use p.gnames and p.vnames  for grid and velocity names
 read filename rhoname, and prop names until prop= "" (up to 20 names)
 output area flowp flowm flow  and for each prop averages */

/* also save to an array  flow,  flowave for each property, for each i-plane 
   aveiflow(noprop+1,iplanes)
 */
void pc_iplaneint(FILE *fpin, FILE *fprint)
{ int *i4d; double *a[4],*x[3],*rho,*u[3];/* need */
  char **gnames,**vnames;  /* p.items  needed */
  
  char *filename,*rhoname,*name[20];  /* input parameters */
  int iname;
  FILE *fp;
  double averages[147]; /* for areaflowint 7*21 */
  double *aveiflow; /* create */
  int irange[4][2]; 
  int iave[4]={0,4,5,6};   /* use from areaflowint */
  char append[4][3]={"a","fp","fm","f"};
  
  int i,j,L,iall;
  double *prop[20];  /* properties */
  
  /* grid and velocities */
  gnames=(char **)need("p.gnames");
  i4d=(int*)need(gnames[0]);
  iall=Prod4(i4d);
  x[0]=(double*)need(gnames[1]); x[1]=x[0]+iall; x[2]=x[1]+iall;
  a[0]=(double*)need(gnames[2]); a[1]=a[0]+i4d[0]; a[2]=a[1]+i4d[1];
  vnames=(char **)need("p.vnames");
  for (j=0;j<3;j++) u[j]=(double*)need(vnames[j]);
  
  /* get input */
  filename=readname(fpin);
  rhoname=readname(fpin);
  rho=(double *)need(rhoname);
  
  printout("normal","iplaneint to file %s, using rho=%s, properties:\n",filename,rhoname);
  iname=0;
  Loop(i,0,20)   /* read and find properties */
  { 
    name[iname]=readname(fpin);
    if (name[iname][0]=='\0')   break;
    if (strcmp(name[iname],"c:")==0) { fseek(fpin,(long)(-3),1);  break; }
    if (arraysize(name[iname])!=iall)
    {
      printout("normal"," %s ignored,",name[iname]);
      continue;
    }
    prop[iname]=(double *)need(name[iname]);
    printout("normal"," %s,",name[iname]);
    iname++;
  }
  printout("normal","\n");
  Loop(i,0,4) {irange[i][0]=0; irange[i][1]=i4d[i]-1; }
  
  aveiflow=(double *)createarray("aveiflow",(iname+1)*i4d[0],'d',0);
  
  fp=safefopen(filename,"w");
  fprintf(fp,"iplaneint \n%d\n i a area flowp flowm flow \n",(iname+1)*4+2);
  for (i=0;i<iname;i++) fprintf(fp," %s%s %s%s %s%s %s%s\n",name[i],append[0],
                                name[i],append[1],name[i],append[2],name[i],append[3]);
  L=0;
  Loop(irange[L][0],0,i4d[L])
  {
    fprintf(fp," %d %lg",irange[L][0],a[L][irange[L][0]]);
    areaflowint(L,irange,i4d,x,rho,u,prop,iname,averages);
    for (i=0;i<=iname;i++) 
    {
      fprintf(fp," %lg %lg %lg %lg\n",averages[7*i+iave[0]],
                                   averages[7*i+iave[1]],averages[7*i+iave[2]],
                                   averages[7*i+iave[3]]);
      aveiflow[i+(iname+1)*irange[L][0]]=averages[7*i+iave[3]];
    }
  }
  fclose(fp);
}