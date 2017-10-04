/* contains c_gridvarmod */
#include "global.h"

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* modify the grid and on-the-points variables together using a b c d, 
 routine will stop if variables are not found or of wrong size */
void c_gridvarmod(FILE *fpin, FILE *fprint)
{   
  double *abcd[4]; int *i4d;   /* need */
  double *xyz;
  char *name[100]; /* names of variables to be modified */
  double *v[100]; /* variables to be modified */
  char cmod,fororo; double aa[2][3];  /* input */
  int iaa[2][3],ifroms,ifrome,nyin;
  int nmod,size,irep;
  int L,L2,L3;
  int iall,ii[4],i,j,k,n,m;
  char cmat[4]={'a','b','c','d'};
  double faa,frac;
  
  abcd[0]=(double *)need("abcd");
  i4d=(int *)need("idim4d");
  xyz=(double *)need("xyz");
  Loop(i,1,4) abcd[i]=abcd[i-1]+i4d[i-1];
  iall=Prod4(i4d);
  nmod=3; v[0]=xyz; v[1]=xyz+iall; v[2]=xyz+2*iall;
  printout("normal","modify: ");
  while(1)  /* read names of variables to be interpolated */
  {
    name[nmod]=readname(fpin);
    if (name[nmod][0]=='\0') break;
    v[nmod]=(double *)need(name[nmod]);
    if (v[nmod]==0) 
    { 
      printout("error c_gridvarmod","%s not found, exit\n",name[nmod]); 
      exitm4d(0); 
    }
    size=arraysize(name[nmod]);
    if (size%iall !=0) 
    {
      printout("error c_gridvarmodl","%s, size %d, not an on-the-points variable, exit\n"
              ,name[nmod],size); 
      exitm4d(0); 
    }
    irep=size/iall;
    printout("normal","%s %d, ",name[nmod],irep);
    k=nmod;
    Loop(i,0,irep)
    { v[k+i]=v[k]+i*iall; nmod++; }
    /*printout("normal"," %d\n",nmod); */
    if (nmod>100) 
    {
      printout("error c_gridvarmod"," %d is too many variables, max is 100\n",nmod); 
      exitm4d(0); 
    }
  }
  printout("normal","\n");
  
  while(1) /* read interp info and do it */
  { 
    cmod=read1charname(fpin);
    L=-1;
    Loop(i,0,3) if (cmod==cmat[i]) L=i;
    if (L==-1)
    {
      printout("normal","end grid modification, cmod %c\n",cmod);
      break;
    }
    L2=(L+1)%3; L3=(L+2)%3;
    
    fororo=read1charname(fpin);
    printout("normal"," move %c",cmod);
    if (fororo=='f') printout("normal"," for region,");
    else printout("normal"," omit region,");
    
    Loop(i,0,3)
    {
      Loop(j,0,2)
      { 
        aa[j][i]=readdouble(fpin);
        iaa[j][i]=findex(aa[j][i],abcd[i],i4d[i],&faa); 
        iaa[j][i]+=faa+.5; 
      }
      if (i==L)
        printout("normal"," %c=%lg index %d dir %lg,",
                cmat[i],aa[0][i],iaa[0][i],aa[1][i]); 
      else
        printout("normal"," %c=%lg,,%lg index=%d,,%d",
                cmat[i],aa[0][i],aa[1][i],iaa[0][i],iaa[1][i]);
    }
    printout("normal","\n");	
    
    frac=aa[1][L];
    if (frac<-1) frac=-1;  else if (frac>1) frac=1;
    /*printout("normal","frac %lg\n",frac); */
    
    Loop(ii[3],0,i4d[3]) Loop(ii[L2],0,i4d[L2]) Loop(ii[L3],0,i4d[L3])
    { 
      nyin=0;
      if (ii[L2]>=iaa[0][L2] && ii[L2]<=iaa[1][L2]
          && ii[L3]>=iaa[0][L3] && ii[L3]<=iaa[1][L3]) nyin=1;
      if (fororo=='f' && nyin==0) continue;
      else if (fororo!='f' && nyin==1) continue;
      ifroms=iaa[0][L]; ifrome=ifroms;
      ii[L]=iaa[0][L];
      m=In4(ii,i4d);
      
      if (frac==1 || frac==-1)  /* find range of identical points */
      {
        for (ii[L]=iaa[0][L]-1; ii[L]>=0; ii[L]--)
        { 
          j=-1;
          n=In4(ii,i4d);
          Loop(k,0,3) if (xyz[m+iall*k]!=xyz[n+iall*k]) {j=0; break; }
          /*  printout("normal","-- i j %d %d\n",i,j); */
          if (j==0) break;
          ifroms--;
        }
        for (ii[L]=iaa[0][L]+1; ii[L]<i4d[L]; ii[L]++)
        {  
          j=1;
          n=In4(ii,i4d);
          Loop(k,0,3) if (xyz[m+iall*k]!=xyz[n+iall*k]) {j=0; break; }
          /*  printout("normal","++ i j %d %d\n",i,j); */
          if (j==0) break;
          ifrome++;
        }
        if (frac<0) ii[L]=max(0,ifroms-1);
        else ii[L]=min(i4d[L]-1,ifrome+1);
        n=In4(ii,i4d);
        /*printout("normal","ifroms ifrome n %d %d %d\n",ifroms,ifrome,n); */
        for (ii[L]=ifroms; ii[L]<=ifrome; ii[L]++)
        {
          m=In4(ii,i4d);
          Loop(k,0,nmod) v[k][m]=v[k][n]; 
        }
      }
      
      else /* fractional move  do just the pt itself */
      { 
        if (frac<0) ii[L]=max(0,ifroms-1);
        else ii[L]=min(i4d[L]-1,ifrome+1);
        n=In4(ii,i4d);
        Loop(k,0,nmod) v[k][m]+=abs(frac)*(v[k][n]-v[k][m]);
      }
    }
  }
}

