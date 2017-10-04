/* contains c_contrhsu */
#include "global.h"

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

/* calc  right hand side continuity, rhsc =  minus the integral  rho Ui Ai,   */
void c_contrhsu(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep,*noindppts;	 	/* needed arrays */ 
  double *u[3],*rho;
  char *cltp;
  int *itrange; /* use if available */
  double *rhs;    /* created or updated array rhsc */
  
  int i4dp[4],ip[4],i,iall,ieq,L,L2,L3,iptc,ig[4],ia[3],sign;
  double area[4][3],ug[3];
  int its,ite;
  int db=0; /* set to 1 for debug prints */
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iall=Prod4(i4d);
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  u[0]=(double *)need("U1");
  u[1]=(double *)need("U2");
  u[2]=(double *)need("U3");
  rho=(double *)need("rho");
  geomcinit();
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  noindppts=(int *)need("noindppts");
  
  rhs=(double *)find("rhsc");
  if (rhs==0) 
  { 
    rhs=(double *)createarray("rhsc",noindppts[0],'d',0);
	 Loop(i,0,noindppts[0]) rhs[i]=0;
  }
  else Loop(i,noindppts[its+1],noindppts[ite+2]) rhs[i]=0;
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) 
  { 
    ieq=wherep[In4(ip,i4dp)];
	 if (db>0)   printout("normal"," ip %d %d %d ieq %d index %d\n",ip[0],ip[1],ip[2],ieq,In4(ip,i4dp));
    if (ieq<0) continue;
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
	 /* for active volumes */
    { 
      rhs[ieq]=0;
      ig[3]=ip[3];
      Loop(L,0,3)   /* each direction */
      { 
        L2=(L+1)%3, L3=(L+2)%3;
        Loop(ia[L],0,2) /* each surface */
        {  
          Loop(i,0,3) ig[i]=ip[i]-1;     /* get areas */
          ig[L]=ig[L]+ia[L];
          geomcavector(area,ig,L);
          Loop(ia[L2],0,2) Loop(ia[L3],0,2) /* indices of 4 corners */
          { 
            Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
            iptc=In4(ig,i4d);
            sign=2*ia[L]-1;
            Loop(i,0,3) ug[i]=u[i][iptc];
            rhs[ieq]-=sign*Dot(ug,area[ia[L2]+2*ia[L3]])*rho[iptc];
            if (db>0)	 Loop(i,0,3) printout("normal"," %d  da ug  %lg %lg\n",i,area[ia[L2]+2*ia[L3]][i],ug[i]);
            if (db>0)	 printout("normal"," ia %d %d %d rhs %lg\n",ia[0],ia[1],ia[2],rhs[ieq]);
          }
        }
      }
    }
  }
}
