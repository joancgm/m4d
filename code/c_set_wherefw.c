/* contains c_set_wherefw */
/* set   noindfwpts = number of independent flow and wall points  
 and wherefw[i] which says where each gridpoint is in this list. 
 Sets wherefw[i]=-1 if it is not a point type which is considered */

#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)

void c_set_wherefw(FILE *fpin, FILE *fprint)
{
  int *match,*i4d,**whoelse;  char *clt; /* needed arrays */
  int *noindfwpts, *wherefw, *whoisfw;   /* create or reset */
  int ninout; char *omit;  /* input parameters */
  int iall,i,j,k,it,nyin,nyptin;
  
  ninout=readint(fpin); 
  nyin=1; if (ninout<0) { nyin=0; ninout=-ninout; }
  omit=(char *)tmalloca(ninout,'c');
  printout("normal","form wherefw ");
  if (nyin==1) printout("normal","including point types ");
  else printout("normal","excluding point types ");
  Loop(i,0,ninout) 
  { 
    omit[i]=read1charname(fpin); 
    printout("normal"," %c",omit[i]); 
  }
  printout("normal","\n");
  
  i4d=(int *)need("idim4d");
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  match=(int *)need("match");
  whoelse=(int **)need("whoelse");
  clt=(char *)need("clt");
  
  noindfwpts=(int *)createarray("noindfwpts",i4d[3]+2,'i',0);
  wherefw=(int *)createarray("wherefw",iall,'i',0);
  Loop(i,0,iall) wherefw[i]=-1;   /* initializa as no equation point  */
  k=0;
  Loop(i,0,i4d[3]+2) noindfwpts[i]=0;
  Loop(i,0,iall)
  { 
    if (nyin==1) 
    { 
      nyptin=0;
      Loop(j,0,ninout) 
      {  
        if (clt[i]==omit[j]) {nyptin=1; break; }
      }
    }
    else        /* exclude */
    { 
      nyptin=1;
      Loop(j,0,ninout) 
      {  
        if (clt[i]==omit[j]) {nyptin=0; break; }
      }
    }
    if (nyptin==1)
    { 
      if (match[i]==i)
      { 
        wherefw[i]=k;
        if (whoelse[i]>0)
        {
          Loop(j,0,whoelse[i][0]) wherefw[whoelse[i][j+1]]=k;
        }
        k++;
        it=i/(i4d[0]*i4d[1]*i4d[2]);
        noindfwpts[it+2]=k;
      }
    }
  }
  
  noindfwpts[0]=k;
  printout("normal"," %d independent points, each t",k);
  Loop(i,1,i4d[3]+2) printout("normal"," %d",noindfwpts[i]); 
  printout("normal","\n");
  
  whoisfw=(int *)createarray("whoisfw",noindfwpts[0],'i',0);
  Loop(i,0,iall) if (wherefw[i]>=0) whoisfw[wherefw[i]]=match[i];
}

