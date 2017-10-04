/* contains c_varinit */
#include "global.h"

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])

/* initialize  variables  */
void c_varinit(FILE *fpin, FILE *fprint)
{ 
  char *name,*clt;  double *abcd;
  double *var,cf,cw,cs;
  int *i4d,i,iall; char iformat;
  char aname[4]={'a','b','c','d'};
  
  i4d=(int *)need("idim4d");
  clt=(char *)need("clt");
  abcd=(double *)need("abcd");
  
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  /* start with on points arrays */
  name=readname(fpin);
  iformat=read1charname(fpin);
  
  if (iformat=='0')    /* read constant */
  { 
    var=(double*)createarray(name,iall,'d',0);
    printout("normal","set %s: format %c",name,iformat);
    cf=readdouble(fpin);
    cw=readdouble(fpin);
    cs=readdouble(fpin);
    printout("normal","  flow %lg, wall %lg, solid %lg\n",cf,cw,cs);
    for (i=0;i<iall;i++) 
    {
      var[i]=cf;
      if (clt[i]=='w') var[i]=cw;
      else if (clt[i]=='s') var[i]=cs;
    }
    return;
  }
  if (iformat=='1') /* prod profile */
  { 
    double *a[4],*aa[4],fac,*v[4],vv,f,vf; int na[4],j,ig[4],k,nynew,nyset;
    
    var=(double*)find(name);
    nynew=0;
    if (var==0)
    {
      var=(double*)createarray(name,iall,'d',0);
      nynew=1;
    }
    printout("normal","set %s: format %c",name,iformat);
    a[0]=abcd;
    Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
    fac=readdouble(fpin);
    printout("normal"," fac=%lg\n",fac);
    Loop(i,0,4) 
    {
      na[i]=readint(fpin);
      aa[i]=(double*)tmalloca(na[i],'d');
      v[i]=(double*)tmalloca(na[i],'d');
    }
    Loop(i,0,4) Loop(j,0,na[i]) aa[i][j]=readdouble(fpin);
    Loop(i,0,4) Loop(j,0,na[i]) v[i][j]=readdouble(fpin);
    Loop(i,0,4) 
    { 
      if (na[i]==1) printout("normal"," %c factor: %lg\n",aname[i],v[i][0]);
		else
		{ 
        printout("normal"," %c ",aname[i]);
        Loop(j,0,na[i]) printout("normal"," %lg",aa[i][j]);
        printout("normal","\n factors:");
        Loop(j,0,na[i]) printout("normal"," %lg",v[i][j]);
        printout("normal","\n");
		}
    }
    Loop(ig[3],0,i4d[3]) Loop3(ig,0,i4d)
    { 
      vv=fac;
      nyset=1;
      Loop(i,0,4) 
      {
        if (na[i]>1) 
        { 
          if (nynew==0 && (a[i][ig[i]]<aa[i][0] || aa[i][ig[i]]>aa[i][na[i]-1])) nyset=0;
          k=findex(a[i][ig[i]],aa[i],na[i],&f);
          vf=v[i][k]+f*(v[i][k+1]-v[i][k]);
          vv*=vf;
        }
        else vv*=v[i][0];
      }
      if (nyset==1) var[In4(ig,i4d)]=vv;
    }
    return;
  }
  if (iformat=='2' || iformat=='3') /* read arrays if 2 fill covered part of array, if 3 fill all */
  { 
    FILE *fp; 
    int i4dold[4],iallold=0,nyi4d=0,nyabcd=0,length,irep,ig[4],ii[2][4],ia[4],mm,m,n,ifw;
    double *aa[4],*a[4],*vin,f[2][4]; char type;
    ifw=iformat;
    a[0]=abcd;
    Loop(i,1,4) a[i]=a[i-1]+i4d[i-1];
    printout("normal"," opening file %s, fill format %c\n",name,iformat);
    fp=safefopen(name,"r");
        
    while(1)
    { 
      name=readname(fp);
      if (name==NULL) break;
      iformat=ifw;
      printout("normal"," %s",name);
      length=readint(fp); 
      type=read1charname(fp);
      printout("normal"," length %d type %c",length,type);
      if (strcmp(name,"idim4d")==0)
      { 
        Loop(i,0,4) i4dold[i]=readint(fp);
        nyi4d=1;
        printout("normal"," array dimensions %d %d %d %d\n",
                i4dold[0],i4dold[1],i4dold[2],i4dold[3]);
        iallold=i4dold[0]*i4dold[1]*i4dold[2]*i4dold[3];
        continue;
      }
      vin=(double *)tmalloca(length,'d');
      Loop(i,0,length) vin[i]=readdouble(fp);
      if (nyi4d==0) 
      { 
        printout("normal", ", can't use do not yet have old grid dimensions, continue\n");
        continue;
      }
      if (strcmp(name,"abcd")==0)
      { 
        aa[0]=vin;
        Loop(i,1,4) aa[i]=aa[i-1]+i4dold[i-1];
        nyabcd=1;
        continue;
      }
      if (nyabcd==0)
      { 
        printout("normal", ", can't use do not yet have old abcd, continue\n");
        continue;
      }
      if (length%iallold !=0)
      {
        printout("normal", ", not an on-the-points array, omit\n");
        continue;
      }
      irep=length/iallold;
      var=(double *)find(name);
      if (var==0) 
      { 
        var=(double*)createarray(name,iall*irep,'d',0);
        if (iformat=='2')
        {
          printout("normal"," created %s, format changed to full fill\n",name);
          iformat='3';
        }
      }
      n=0;
      Loop(ig[3],0,i4d[3]) Loop3(ig,0,i4d)
      {  
        if (iformat=='2')  /* omit outside of range with this option */
        { 
          if (a[0][ig[0]]<aa[0][0] || a[0][ig[0]]>aa[0][i4dold[0]-1] ) continue;
          if (a[1][ig[1]]<aa[1][0] || a[1][ig[1]]>aa[1][i4dold[1]-1] ) continue;
          if (a[2][ig[2]]<aa[2][0] || a[2][ig[2]]>aa[2][i4dold[2]-1] ) continue;
          if (a[3][ig[3]]<aa[3][0] || a[3][ig[3]]>aa[3][i4dold[3]-1] ) continue;
        }
        Loop(i,0,4) 
        {
          if (i4dold[i]==1) {ii[0][i]=0; f[1][i]=0; }
          else ii[0][i]=findex(a[i][ig[i]],aa[i],i4dold[i],&f[1][i]);
          ii[1][i]=min(ii[0][i]+1,i4dold[i]-1);
          f[0][i]=1.-f[1][i];
        }
        mm=In4(ig,i4d);
        Loop(i,0,irep) var[mm+i*iall]=0;
        Loop(ia[0],0,2) Loop(ia[1],0,2) Loop(ia[2],0,2) Loop(ia[3],0,2)
        {
          m=ii[ia[0]][0]+i4dold[0]*(ii[ia[1]][1]+i4dold[1]*
                                    (ii[ia[2]][2]+i4dold[1]*ii[ia[3]][3]));
          Loop(i,0,irep)																  
          var[mm+i*iall] += f[ia[0]][0]*f[ia[1]][1]*f[ia[2]][2]*f[ia[3]][3]*vin[m+i*iallold];
        }
        n+=irep;
      }
      printout("normal"," set %d of %d values\n",n,iall*irep);
    }
  }
}
