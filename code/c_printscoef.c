/* contains c_printscoef */
/* print selected coefficients from either a on the grid or pc grid coefficient matrix */ 
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_printscoef(FILE *fpin, FILE *fprint)
{ 
  int *i4d, *match=0, *where, *cn, **ci;  /* needed (maybe) arrays */
  double **cc;
  
  char *namen, *namei, *namec;    /* input parameters */
  int nypc,iexi,jrep,is[4],ie[4];
  
  int i,j,ideq[4],idci[4],ii[4],ipt,ieq,nyprint,ialt[4],nyin;
  double csum;
  
  namen=readname(fpin);
  namei=readname(fpin);
  namec=readname(fpin);
  nypc=readint(fpin);
  iexi=readint(fpin);
  jrep=readint(fpin);
  Loop(i,0,4) {is[i]=readint(fpin); ie[i]=readint(fpin); }
  printout("print"," %s %s %s, nypc %d, iexi %d, jrep %d, range i %d %d j %d %d k %d %d t %d %d\n",
          namen,namei,namec,nypc,iexi,jrep,is[0],ie[0],is[1],ie[1],is[2],ie[2],is[3],ie[3]);
  
  cn=(int *)find(namen);
  ci=(int **)find(namei);
  cc=(double **)find(namec);
  if (cn==0 || ci==0 || cc==0)
  { 
    printout("print","needed array not found %s %s or  %s\n",namen,namei,namec);
    return; 
  }
  i4d=(int *)need("idim4d");
  Loop(i,0,4) { ideq[i]=i4d[i]; idci[i]=i4d[i]; }
  if (nypc==1) 
  {
    Loop(i,0,3) ideq[i]++;
    where=(int *)need("wherep");
  }
  else
  {
    where=(int *)need("wherefw");
    match=(int *)need("match");
  }
  if (iexi==1) Loop(i,0,3) idci[i]++;
  Loop(i,0,4) ie[i]=min(ie[i],ideq[i]);
  
  Loop(ii[3],is[3]-1,ie[3]) Loop(ii[2],is[2]-1,ie[2]) 
  Loop(ii[1],is[1]-1,ie[1]) Loop(ii[0],is[0]-1,ie[0])
  { 
    nyprint=1;
    ipt=In4(ii,ideq);
    ieq=where[ipt];
    if (ieq<0) nyprint=0;
    if (nypc==0 && match[ipt]!=ipt)
    { 
      iexpand(match[ipt],ideq,ialt);
      nyin=1;
      Loop(i,0,4)
      { 
        if (ialt[i]<is[i]-1) nyin=0;
        if (ialt[i]>ie[i]-1) nyin=0;
      }
      if (nyin==0) nyprint=1;
      ieq=where[match[ipt]];
      if (ieq<0) nyprint=0;
    }
    if (nyprint==1)
    { 
      printout("print","at %d %d %d %d, ieq %d, no of coefs %d",
              ii[0],ii[1],ii[2],ii[3],ieq,cn[ieq]);
      if (cn[ieq]>0)
      { 
        printout("print",": sums");
        Loop(j,0,jrep)
        {
          csum=0;
          Loop(i,0,cn[ieq]) csum+=cc[ieq][i+cn[ieq]*j];
          printout("print"," %lg",csum);
        }
        printout("print",", pt, expanded, coefs\n");
        Loop(i,0,cn[ieq])
        { 
          printout("print","           %d",ci[ieq][i]);
          if (iexi==0 || iexi==1) 
          {
            iexpand(ci[ieq][i],idci,ialt);
            if (iexi==nypc) printout("print"," %d %d %d %d",
                                    ialt[0]-ii[0],ialt[1]-ii[1],ialt[2]-ii[2],ialt[3]-ii[3]);
            else  printout("print"," %d %d %d %d",
                          ialt[0],ialt[1],ialt[2],ialt[3]);
          }
          Loop(j,0,jrep) printout("print","  %lg",cc[ieq][i+cn[ieq]*j]);
          printout("print","\n");
        }
      }
      else {printout("print","\n"); }
    }
  }
}
