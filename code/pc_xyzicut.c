/* contains pc_xyzicut */
#include "global.h"
/* cut the grid in the i-grid direction
 to form planes of constant ax+by+cz=s[noplanes] 
 choose append name (e.g. cut) then
 create idim4dcut xyzcut abcdcut cltcut  and varcut for all on the points variables listed */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
/* ---------------------------- */
void pc_xyzicut(FILE *fpin, FILE *fprint)
{
  int *i4d; double *a[4],*x[3]; char *clt, *csym; /* need */
  char *append; int nopl; double axyz[3],*s,ss,se,ds;
  double *v;
  double bs[2],be[2]; char *name; /* input parameters */
  int *i4dnew; double *anew[4], *xnew[3]; char *cltnew=0; double *var; /* create */
  
  int iall,i,ireps[2]={0,0},irepe[2]={0,0},jks[2],jke[2],n,iallnew,nn;
  int ir,j,iinew[4],ii[4],nyok;
  int *inew,*jknew[3],*irnew[3]; double *fnew, *si;  /* temporary arrays */
  char namg[4][20]={"idim4d","abcd","xyz","clt"},namenew[50];
  double dxr[3][3]={{0,0,0},{0,0,0},{0,0,0}},xr[3],f;
  char cltp;
  Array *array;
  
  /* get needed arrays */
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  a[0]=(double *)need("abcd");
  Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[1]+iall;
  clt=(char *)need("clt");
  csym=(char *)need("csym");
  /* read input and set parameters */
  append=readname(fpin);
  nopl=readint(fpin);
  
  Loop(i,0,3) axyz[i]=readdouble(fpin);
  printout("normal"," append %s, nopl=%d, %lgx+%lgy+%lgz = s, \n s=",
          append,nopl,axyz[0],axyz[1],axyz[2]);
  if (nopl>0) 
  {
    s=(double *)tmalloca(nopl,'d');
    Loop(i,0,nopl) 
    {
      s[i]=readdouble(fpin);
      printout("normal"," %lg",s[i]);
    }
  }
  else
  {
    ss=readdouble(fpin);
    se=readdouble(fpin);
    ds=readdouble(fpin);
    nopl=(se-ss)/ds+.1;
    s=(double *)tmalloca(nopl,'d');
    Loop(i,0,nopl) s[i]=ss+i*ds;
    printout("normal"," %lg to %lg by %lg",s[0],s[nopl-1],ds);
  }
  printout("normal","\n");
  Loop(i,0,2)
  {
    bs[i]=readdouble(fpin);
    be[i]=readdouble(fpin);
    if (csym[8+2*i]=='r')
    {
      for (n=0;n<10;n++)
      {
        if (bs[i]>=a[i+1][0]) break;
        ireps[i]--;
        bs[i]+=a[i+1][i4d[i+1]-1]-a[i+1][0];
      }
      for (n=0;n<10;n++)
      {
        if (be[i]<=a[i+1][i4d[i+1]-1]) break;
        irepe[i]++;
        be[i]-=a[i+1][i4d[i+1]-1]-a[i+1][0];
      }
    }
    else 
    {
      bs[i]=max(a[i+1][0],bs[i]);
      be[i]=min(a[i+1][i4d[i+1]-1],be[i]);
    }  
    for (n=1;n<i4d[i+1];n++) if (bs[i]<a[i+1][n]) break;
    jks[i]=n-1;
    for (n=i4d[i+1]-2;n>=0;n--) if (be[i]>a[i+1][n]) break;
    jke[i]=n+1;
    printout("normal","dim %d: bc from %lg to %lg,jk from %d irep %d, to %d irep %d\n",
            i+1,bs[i],be[i],jks[i],ireps[i],jke[i],irepe[i]);
  }
  /* end of first part of input, prop names read later */
  /* set up new grid */
  Loop(i,0,4) strncat(namg[i],append,strlen(append)+1);
  i4dnew=(int *)createarray(namg[0],4,'i',1);
  i4dnew[0]=nopl;
  Loop(i,1,3) i4dnew[i]=jke[i-1]-jks[i-1]+1+(i4d[i]-1)*(irepe[i-1]-ireps[i-1]);
  i4dnew[3]=i4d[3];
  iallnew=Prod4(i4dnew);
  printout("normal","%s = %d %d %d %d\n",namg[0],i4dnew[0],i4dnew[1],i4dnew[2],i4dnew[3]);
  printout("normal"," creating %s %s %s\n",namg[1],namg[2],namg[3]);
  /* set abcdnew */
  anew[0]=(double *)createarray(namg[1],i4dnew[0]+i4dnew[1]+i4dnew[2]+i4dnew[3],'d',1);
  Loop(i,1,4) anew[i]=anew[i-1]+i4dnew[i-1];
  Loop(i,0,nopl) anew[0][i]=s[i];
  jknew[1]=(int *)tmalloca(i4dnew[1],'i');
  jknew[2]=(int *)tmalloca(i4dnew[2],'i');
  irnew[1]=(int *)tmalloca(i4dnew[1],'i');
  irnew[2]=(int *)tmalloca(i4dnew[2],'i');
  
  Loop(n,0,2) 
  {
    ir=ireps[n];
    j=jks[n];
    Loop(i,0,i4dnew[n+1])
    {
      anew[n+1][i]=a[n+1][j]+ir*(a[n+1][i4d[n+1]-1]-a[n+1][0]);
      jknew[n+1][i]=j;
      irnew[n+1][i]=ir;
      j++;
      if (j>i4d[n+1]-1) {ir++; j=1;}
    }
  }
  Loop(i,0,i4d[3]) anew[3][i]=a[3][i];
  /* set xyznew cltnew(as old but 's' in not found in i-range
   and temp arrays inew[i,j,k] fnew[i,j,k] locating points in new grid to interpolate properties to new grid */
  xnew[0]=(double *)createarray(namg[2],iallnew*3,'d',1);
  xnew[1]=xnew[0]+iallnew; xnew[2]=xnew[1]+iallnew;
  cltnew=(char *)createarray(namg[3],iallnew,'c',1);
  fnew=(double *)tmalloca(iallnew,'d');
  inew=(int *)tmalloca(iallnew,'i');
  si=(double *)tmalloca(i4d[0],'d');
  Loop(iinew[3],0,i4d[3])
  {
    Loop(i,0,3)
    {
      dxr[1][i]=x[i][i4d[0]*(i4d[1]-1)]-x[i][0];
      dxr[2][i]=x[i][i4d[0]*i4d[1]*(i4d[2]-1)]-x[i][0];
    }
    Loop(iinew[1],0,i4dnew[1]) Loop(iinew[2],0,i4dnew[2])
    { 
      ii[1]=jknew[1][iinew[1]];
      ii[2]=jknew[2][iinew[2]];
      ii[3]=iinew[3];
      Loop(ii[0],0,i4d[0])
      { 
        n=In4(ii,i4d);
        Loop(i,0,3) xr[i]=x[i][n]+irnew[1][iinew[1]]*dxr[1][i]+irnew[2][iinew[2]]*dxr[2][i];
        si[ii[0]]=axyz[0]*xr[0]+axyz[1]*xr[1]+axyz[2]*xr[2];
      }
     /*  printout("normal"," si-range %lg %lg\n",si[0],si[i4d[0]-1]); */
      Loop(iinew[0],0,i4dnew[0])
      {
        nn=In4(iinew,i4dnew);
       /*  printout("normal"," iinew %d %d %d %d, nn %d \n",iinew[0],iinew[1],iinew[2],iinew[3],nn); */
        cltnew[nn]='s'; 
        fnew[nn]=0; 
        ii[0]=0; 
        nyok=0; 
        if (si[0]==si[i4d[0]-1]) nyok=1;
        else if ((s[iinew[0]]-si[0])/(si[i4d[0]-1]-si[0]) < 0) nyok=1;
        else if ((s[iinew[0]]-si[0])/(si[i4d[0]-1]-si[0]) > 1) 
        { cltnew[nn]='s'; fnew[nn]=1; ii[0]=i4d[0]-2; nyok=1;}
        else
        { 
          Loop(ii[0],0,i4d[0]-1)
          {
            if (si[ii[0]+1]==si[ii[0]]) continue;
            f=(s[iinew[0]]-si[ii[0]])/(si[ii[0]+1]-si[ii[0]]);
            if (f>-1.e-6 && f<1.+1e-6)
            { 
              nyok=1;
              fnew[nn]=max(0,min(f,1));
              cltnew[nn]=clt[In4(ii,i4d)];
              if (fnew[nn]==0) break;
              cltp=clt[In4(ii,i4d)+1];
              if (fnew[nn]==1) { cltnew[nn]=cltp; break;}
              if (cltp=='s') break;
              if (cltnew[nn]=='s') { cltnew[nn]=cltp; break;}
              if (cltp=='w') break;
              if (cltnew[nn]=='w') { cltnew[nn]=cltp; break;}  
              break;
            }
          }
        }
        /* printout("normal"," nyok %d ii %d %d %d %d fnew %lg cltnew %c\n",
               nyok,ii[0],ii[1],ii[2],ii[3],fnew[nn],cltnew[nn]);  */
        n=In4(ii,i4d);
        inew[nn]=n;
        Loop(i,0,3) xnew[i][nn]=
        irnew[1][iinew[1]]*dxr[1][i]+irnew[2][iinew[2]]*dxr[2][i]
        +(1.-fnew[nn])*x[i][n]+fnew[nn]*x[i][n+1];
        /* printout("normal","xnew %lg %lg %lg\n",xnew[0][nn],xnew[1][nn],xnew[2][nn]); */
      }
    }
  }
  /* read properties to be converted */
  printout("normal"," converting arrays:\n");
  n=100;
  while (n)
  { 
    name=readname(fpin);
    if (name[0]=='\0') break;
    if (strcmp(name,"c:")==0) break;
    n--;
    array=findarray(name);
    if (array==0)  {printout("warning pc_xyzicut"," array %s not found,\n",name); continue; }
  
    if (array->type!='d' || (array->size)%iall !=0) 
    { printout("warning pc_xyzicut"," not on points type d, %s ignored\n",name); continue; }
    printout("normal"," %s",name);
    j=(array->size)/iall;
    v=(double *)array->pointer;
    strcpy(namenew,name);
    strncat(namenew,append,strlen(append)+1);
    var=(double *)createarray(namenew,j*iallnew,'d',1);
    printout("normal","  %s\n",namenew);
    Loop(iinew[3],0,i4d[3])
    Loop3(iinew,0,i4dnew)
    {
      nn=In4(iinew,i4dnew);
     /*  printout("normal","%d %d %d %d nn %d\n",iinew[0],iinew[1],iinew[2],iinew[3],nn);
      printout("normal","f %lg inew %d v %lg %lg\n",fnew[nn],inew[nn]);  */
      Loop(i,0,j) var[nn+i*iallnew]=(1-fnew[nn])*v[inew[nn]+i*iall]+fnew[nn]*v[inew[nn]+1+i*iall];
    }
  }
}