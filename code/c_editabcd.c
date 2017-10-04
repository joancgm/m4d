/* contains c_editabcd */
#include "global.h"
/* edit on the points arrays (but with/without jrep) using abcd */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_editabcd(FILE *fpin, FILE *fprint)
{ 
  int *i4d; double *a[4]; /* needed arrays */
  char *name,type,action;  /* input parameters */
  int jrep,js,je;
  double ar[4][2];
  
  Array *af; 
  int i,j,iall,is[4],ie[4],ii[4],ipt,k;
  int *mf=0,mv=0; double *df=0,dv=0,f; char *cf=0,cv=' ';
  char nabcd[4]={'a','b','c','d'};
  char ijkn[4]={'i','j','k','n'};
  
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double *) need("abcd");
  Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
  
  name=readname(fpin);     /* array info  */
  type=read1charname(fpin);  
  jrep=readint(fpin);
  if (type !='d' && type !='c' && type !='i')
  { 
    printout("error c_editabcd","error, can't edit array of type %c\n",type); 
    exitm4d(0); 
  }
  
  af=findarray(name);
  if (af==0)   /* create and clear array */
  { 
    if (type=='c') 
    {
      cf=(char *)createarray(name,iall*jrep,'c',1); 
      Loop(i,0,iall*jrep) cf[i]=(char)0; 
    }
    else if (type=='i') 
    { 
      mf=(int *)createarray(name,iall*jrep,'i',1); 
      Loop(i,0,iall*jrep) mf[i]=0; 
    }
    else if (type=='d') 
    { 
      df=(double *)createarray(name,iall*jrep,'d',1); 
      Loop(i,0,iall*jrep) df[i]=0; 
    }
    printout("normal","array %s created and cleared, type=%c, jrep=%d, size %d\n",name,type,jrep,iall*jrep);
  }
  else   /* check */
  { 
    if (type != af->type) 
    {
      printout("error c_editabcd","error can't edit %s, type error %c %c\n",name,type,af->type); 
      exitm4d(0); 
    }
    if (iall*jrep != af->size)
    {
      printout("error c_editabcd","error can't edit %s, size error %d %d\n",name,iall*jrep,af->size); 
      exitm4d(0); 
    }
    if (type=='c') cf=(char *)af->pointer;
    else if (type=='i') mf=(int *)af->pointer;
    else if (type=='d') df=(double *)af->pointer;
    printout("normal","array %s found,  type=%c, jrep=%d, size %d\n",name,type,jrep,iall*jrep);
  }
  
  /* ok have array, now specify what to do, repeat until action is not valid, not  s-set, a-add, m-multiply */
  while(1)
  {  
    action=read1charname(fpin);
    printout("normal"," %c",action); 
    if (type=='c' && action !='s') {printout("normal","  end edit\n"); return;}
    if (action!='s' && action !='a' && action !='m') {printout("normal"," end edit\n"); return;}
    if (type=='c') {cv=read1charname(fpin); printout("normal","  %c   ",cv); }
    else if (type=='d') {dv=readdouble(fpin); printout("normal"," %lg    ",dv);}
    else if (type=='i') {mv=readint(fpin); printout("normal"," %d   ",mv); }
    Loop(i,0,4)  /* read abcd */
    { 
      ar[i][0]=readdouble(fpin);
      ar[i][1]=readdouble(fpin);
      if (i4d[i]==1) { is[i]=0; ie[i]=1; }
      else 
      { 
        is[i]=findex(ar[i][0],a[i],i4d[i],&f);
        if (f>.9) is[i]++;
        ie[i]=1+findex(ar[i][1],a[i],i4d[i],&f);
        if (f>.9) ie[i]++;
      }
      printout("normal","   %c %lg %lg (%c %d %d)",nabcd[i],ar[i][0],ar[i][1],ijkn[i],is[i],ie[i]-1);
    }
    if (jrep==1) {js=0; je=1;}
    else   { 
      js=readint(fpin);
      je=max(readint(fpin),js+1);
      printout("normal","jrep %d %d",js,je);
    }
    printout("normal","\n");
    Loop(ii[0],is[0],ie[0]) Loop(ii[1],is[1],ie[1]) Loop(ii[2],is[2],ie[2]) Loop(ii[3],is[3],ie[3]) 
    { 
      ipt=In4(ii,i4d);
      Loop(j,js,je) 
      { 
        k=ipt+j*iall;
        if (type=='c') cf[k]=cv;
        else if (type=='d')
        { 
          if (action=='s') df[k]=dv;
          else if (action=='a') df[k]+=dv;
          else if (action=='m') df[k]*=dv;
        }
        else if (type=='i')
        { 
          if (action=='s') mf[k]=mv;
          else if (action=='a') mf[k]+=mv;
          else if (action=='m') mf[k]*=mv;
        }
      }
    }
  }
}
