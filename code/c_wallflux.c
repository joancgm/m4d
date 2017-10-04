/* contains c_wallflux */
/* add to rhsin the specified wall flux per unit area times the area */
/* for this first version fluxda is an array of size 1 */
 
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
/* iexpand (k,idim,i) int k,idim[4],i[4]; 
 Expand index k into its 4 components i. based on dimensions idim
 */
 
/*----------------------------------------------*/
void c_wallflux(FILE *fpin, FILE *fprint)
{ 
  int *i4d, *noindwpts, *whoisw, **whoelse, *wherefw;  /* need */
  char *clt; 
  double *xyzd; 
  char *namerhsin, *namefluxda;   /* read */
  double *rhsin, *fluxda;
  
  int i,j,k,iwi,iw,i4dd[4],ialld,ii[4],id[4],iip[4],klim;
  int mg[3][3][3],md[3][3][3],m,ia[3],L,L2,L3,i2,i3,j2,j3;
  double dy[2][3],flux,area[3];
  char nywall;
    
  namerhsin=readname(fpin);
  namefluxda=readname(fpin);
  rhsin=(double *)need(namerhsin);
  fluxda=(double *)need(namefluxda);
  printout("normal","adding wall flux to %s, flux/area= %s\n",namerhsin,namefluxda);
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dd[i]=i4d[i]*2-1; i4dd[3]=i4d[3];
  ialld=Prod4(i4dd);
  noindwpts=(int *)need("noindwpts");
  whoisw=(int *)need("whoisw");
  whoelse=(int **)need("whoelse");
  wherefw=(int *)need("wherefw");
  
  clt=(char *)need("clt");
  xyzd=(double *)need("xyzdouble");
  
  Loop(j,0,noindwpts[0])    /* each independent wall point */
  {
    iwi=whoisw[j];
    klim=1;
    if (whoelse[iwi]>0) klim=whoelse[iwi][0]+1;
    Loop(k,0,klim)  /* each wall point iw */
    {
      if (k==klim-1) iw=iwi;
      else iw=whoelse[iwi][k+1];
      iexpand(iw,i4d,ii);         /* ijk for wall point */
      Loop(i,0,3) id[i]=2*ii[i]; id[3]=ii[3];  /* ijk for double point */
      /* determine indices of surrounding points */
      iip[3]=ii[3];
      Loop(ia[0],0,3)  Loop(ia[1],0,3)  Loop(ia[2],0,3) 
          {
            for (i=0;i<3;i++) iip[i]=max(0,min(i4d[i]-1,ii[i]+ia[i]-1));
            mg[ia[0]][ia[1]][ia[2]]=In4(iip,i4d);
            for (i=0;i<3;i++) iip[i]=max(0,min(i4dd[i]-1,id[i]+ia[i]-1));
            md[ia[0]][ia[1]][ia[2]]=In4(iip,i4dd);
          }
      /* look at each possible surface direction */
      Loop(L,0,3)
      {
        L2=(L+1)%3; L3=(L+2)%3;
        ia[L]=1;
        Loop(i2,0,2) Loop(i3,0,2) /* each of 4 surfaces */
        {
           nywall='y';
          Loop(j2,0,2) Loop(j3,0,2) /* each of 4 points */
          {
            ia[L2]=i2+j2; ia[L3]=i3+j3;
            m=mg[ia[0]][ia[1]][ia[2]];
            if (clt[m] != 'w') nywall='n';
          }
          if (nywall=='n') continue;
          /* have wall section */
          Loop(i,0,3) { dy[0][i]=0; dy[1][i]=0; }
          Loop(j2,0,2) Loop(j3,0,2) /* each of 4 points */
          {
            ia[L2]=i2+j2; ia[L3]=i3+j3;
            m=md[ia[0]][ia[1]][ia[2]];
            for (i=0;i<3;i++)
            {
              dy[0][i]+=xyzd[m+ialld*i]*(2*j2-1);
              dy[1][i]+=xyzd[m+ialld*i]*(2*j3-1);
            }
          }
          Cross(dy[0],dy[1],area);
          flux=.5*sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
          rhsin[wherefw[iwi]]+=flux*fluxda[0];
        }
      }
    }
  }
}