/* contains c_coefzero */
#include "global.h"
/* set coefs to zero, or center point only  for specified point types  */

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_coefzero(FILE *fpin, FILE *fprint)
{ 
  int *noindfwpts,*cn, **ci,*whoisfw,*nocoefs;   char *clt;/* needed arrays */
  double **cc;  /* modify */
  int *itrange; double *roundoff; /* use if available */
  int jzero,jrep,nycenter1,nexclude;  /* input parameters */
  char *ncn,*nci,*ncc;
  char exclude[20]; 
  
  int i,k,ipt,n;
  double tol;
  int neqs,neqe;
  
  ncn=readname(fpin);
  nci=readname(fpin);
  ncc=readname(fpin);
  jrep=readint(fpin);
  jzero=readint(fpin);
  nycenter1=readint(fpin);
  if (strcmp(ncn,"coef_n")==0)  /* use nocoefs instead of jrep */
  {
    nocoefs=(int *)need("nocoefs");
    jrep=nocoefs[0];
    if (jzero>nocoefs[0]-1)
    { 
      printout("error c_coefzero"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
      exitm4d(0);
    }
  }
  printout("normal"," coefn,i,c %s %s %s jrep %d jzero %d nycenter1 %d\n",
          ncn,nci,ncc,jrep,jzero,nycenter1);
  nexclude=readint(fpin);
  printout("normal"," exclude %d pt types: ",nexclude);
  if (nexclude>20) nexclude=20;
  Loop(i,0,nexclude)
  {
    exclude[i]=read1charname(fpin);
    printout("normal"," %c",exclude[i]);
  }
  printout("normal","\n");
  
  noindfwpts=(int *)need("noindfwpts");/* get needed arrays */
  whoisfw=(int *)need("whoisfw");
  clt=(char *)need("clt");
  cn=(int *)need(ncn);
  ci=(int **)need(nci);
  cc=(double **)need(ncc);
  roundoff=(double *)find("roundoff");
  tol=1.e-8; if (roundoff>0) tol=roundoff[0];
  itrange=(int *)find("itrange");
  if (itrange==0) {neqs=0; neqe=noindfwpts[0]; }
  else {neqs=noindfwpts[itrange[0]+1]; neqe=noindfwpts[itrange[1]+2]; }
  
  Loop(i,neqs,neqe)
  if (cn[i]>0)
  { 
    ipt=whoisfw[i];
    Loop(k,0,nexclude)
    if (clt[ipt]==exclude[k])
    { 
      if (jrep==1) { free(cc[i]); free(ci[i]); cn[i]=0; }
      else
      { 
        Loop(n,0,cn[i]) cc[i][n+jzero*cn[i]]=0;
        if (nycenter1>0)
        {
          Loop(n,0,cn[i]) 
          if (ipt==ci[i][n]) 
          {
            cc[i][n+jzero*cn[i]]=1; 
            break;  
          }
        }
        else coefcombine(i,cn,ci,cc,jrep,tol); 
      }
    }
  }
}


