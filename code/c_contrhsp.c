/* contains c_contrhsp */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

/*  add permanent p contribution to rhsc */

void c_contrhsp(FILE *fpin, FILE *fprint)
{ 
  int *noindppts,*cpc_n,**cpc_i;      /* needed arrays */
  double *pp, **cpc_c;
  int *itrange;
  double *rhsc; /* update */
  int ieq,i,n;
  int neqs,neqe;
  
  noindppts=(int *)need("noindppts");
  cpc_n=(int *)need("cpc_n");
  cpc_i=(int **)need("cpc_i");
  pp=(double *)need("pp");
  rhsc=(double *)need("rhsc");
  cpc_c=(double **)need("cpc_c");
  itrange=(int *)find("itrange");
  if (itrange==0) { neqs=0; neqe=noindppts[0]; }
  else { neqs=noindppts[itrange[0]+1]; neqe=noindppts[itrange[1]+2]; }
  
  Loop(ieq,neqs,neqe)
  { 
    n=cpc_n[ieq];
	 if (n>0)
      Loop(i,0,n) 
      rhsc[ieq]-=pp[cpc_i[ieq][i]]*cpc_c[ieq][i+n];
  }
}
