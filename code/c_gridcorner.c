/* contains c_gridcorner */
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* fix grid at corners using abc values and store array ijkcorner */
void c_gridcorner(FILE *fpin, FILE *fprint)
{  
  int *i4d; double *a[3]; /* needed arrays */
  double *x[3];  /* modify */
  int *ijkcorner; /* create */
  
  char *cname;	/* input */
  double aa[2][3];
  
  int iaa[2][3],is[3],im[3]; 
  double faa[2][3],xave[3];
  int i,j,k,iall,ii[4],ipt[4],n,npt,L;
  char nabc[4]="abc";
  int ncorn,icorn[1000];
  
  i4d=(int *)need("idim4d");
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall;
  x[2]=x[1]+iall;
  a[0]=(double *)need("abcd");
  a[1]=a[0]+i4d[0];
  a[2]=a[1]+i4d[1];
  
  ncorn=0;
  while(1)
  { 
    cname=readname(fpin);
    printout("normal"," %s",cname);
    if (cname[0]=='\0' || cname[0]=='e')
    { printout("normal","end gridcorner\n"); break; }
    if (strcmp(cname,"c:")==0) 
    { fseek(fpin,(long)(-3),1);printout("normal","end gridcorner\n");  break; }
    Loop(i,0,3) 
    { 
      printout("normal",", %c",nabc[i]);
      Loop(j,0,2) 
      { 
        aa[j][i]=readdouble(fpin);
        printout("normal"," %lg",aa[j][i]); 
        iaa[j][i]=findex(aa[j][i],a[i],i4d[i],&faa[j][i]); 
        iaa[j][i]+=faa[j][i]+.001; 
      }
    }
    printout("normal",", i= %d %d, j=%d %d, k=%d %d",
            iaa[0][0],iaa[1][0],iaa[0][1],iaa[1][1],iaa[0][2],iaa[1][2]);
    Loop(i,0,3) {is[i]=1; if (iaa[1][i]<iaa[0][i]) is[i]=-1;}
    Loop(i,0,3) im[i]=is[i]*(iaa[1][i]-iaa[0][i])+1;
    printout("normal",", im %d %d %d\n",im[0],im[1],im[2]);
    
    /* store info */
    if (strcmp(cname,"ab")==0) icorn[ncorn]=1;
    ncorn++;
    Loop(i,0,3) {icorn[ncorn+i]=iaa[0][i]; icorn[ncorn+4+i]=iaa[1][i]; }
    icorn[ncorn+3]=0; icorn[ncorn+7]=i4d[3];
    ncorn+=8;
    if (ncorn>999) 
    {
      printout("error c_gridcorner","too many corners %d, redimension\n",ncorn/9); 
      exitm4d(0); 
    }
    
    if (strcmp(cname,"ab")==0)
    { 	
      if (im[0]!=im[1]) 
      { 
        printout("error c_gridcorner","can't do unequal di dj on corner, exit\n"); 
        exitm4d(0); 
      }
      Loop(i,0,im[0])
      { 
        ii[0]=iaa[0][0]; ii[1]=iaa[0][1]+is[1]*i;
        Loop(k,0,im[2])
        Loop(ii[3],0,i4d[3])
        {  
          ii[2]=iaa[0][2]+is[2]*k;
          ipt[0]=iaa[0][0]+is[0]*i; 
          ipt[1]=iaa[0][1];
          ipt[2]=ii[2]; 
          ipt[3]=ii[3];
          n=In4(ii,i4d); npt=In4(ipt,i4d);
          /*printout("normal"," ii %d %d %d %d n %d xyz %lg %lg %lg\n",
           ii[0],ii[1],ii[2],ii[3],n,x[0][n],x[1][n],x[2][n]);
           printout("normal"," ipt %d %d %d %d npt %d xyz %lg %lg %lg\n",
           ipt[0],ipt[1],ipt[2],ipt[3],npt,x[0][npt],x[1][npt],x[2][npt]);*/
          Loop(L,0,3) 
          { 
            xave[L]=.5*(x[L][n]+x[L][npt]);
            x[L][n]=xave[L]; x[L][npt]=xave[L];
          }
          /*printout("normal"," ii %d %d %d %d n %d xyz %lg %lg %lg\n",
           ii[0],ii[1],ii[2],ii[3],n,x[0][n],x[1][n],x[2][n]);
           printout("normal"," ipt %d %d %d %d npt %d xyz %lg %lg %lg\n",
           ipt[0],ipt[1],ipt[2],ipt[3],npt,x[0][npt],x[1][npt],x[2][npt]);*/
          Loop(j,1,i+1)
          { 
            ipt[0]=iaa[0][0]+is[0]*j;
            ipt[1]=iaa[0][1]+is[1]*i;
            npt=In4(ipt,i4d);
            Loop(L,0,3) x[L][npt]=xave[L];
            /*printout("normal"," ipt %d %d %d %d npt %d xyz %lg %lg %lg\n",
             ipt[0],ipt[1],ipt[2],ipt[3],npt,x[0][npt],x[1][npt],x[2][npt]);*/						
            ipt[0]=iaa[0][0]+is[0]*i;
            ipt[1]=iaa[0][1]+is[1]*j;
            npt=In4(ipt,i4d);
            Loop(L,0,3) x[L][npt]=xave[L];
            /*printout("normal"," ipt %d %d %d %d npt %d xyz %lg %lg %lg\n",
             ipt[0],ipt[1],ipt[2],ipt[3],npt,x[0][npt],x[1][npt],x[2][npt]);*/
          }
        }
      }
    }
    else
    { 
      printout("error c_gridcorner","%s not yet implemented, exit\n",cname);
      exitm4d(0);
    }
  }
  /* perm store corner info */
  icorn[ncorn]=0;
  ncorn++;
  ijkcorner=(int *)createarray("ijkcorner",ncorn,'i',0);
  Loop(i,0,ncorn) ijkcorner[i]=icorn[i];
}
