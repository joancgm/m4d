/* contains c_coefvisc c_coefviscstep coefvisc */
#include "global.h"
/* calculate coefficient for viscous term based on grad at 8 corners
 c_coefvisc interpolates grad to center of surface
 c_coefviscstep uses best combination of two points either side of surface
 for either isotropic or tensor viscosities  */

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 

void coefvisc(int nylin, char *namevisc, int jvisc, FILE *fprint);
/* ---------------------------- */
void c_coefvisc(FILE *fpin, FILE *fprint)
{ int nylin=1;
  char *namevisc; int jvisc;   /* input parameter */
  
  namevisc=readname(fpin);
  jvisc=readint(fpin);
  printout("normal"," viscosity: %s, at coef %d\n",namevisc,jvisc);
  
  coefvisc(nylin,namevisc,jvisc,fprint);
}
/* ---------------------------- */
void c_coefviscstep(FILE *fpin, FILE *fprint)
{ int nylin=0;
  char *namevisc; int jvisc;   /* input parameter */
  
  namevisc=readname(fpin);
  jvisc=readint(fpin);
  printout("normal"," viscosity: %s, at coef %d\n",namevisc,jvisc);
  
  coefvisc(nylin,namevisc,jvisc,fprint);
}
/* ---------------------------- */
void coefvisc(int nylin, char *namevisc, int jvisc, FILE *fprint)
{
  int *i4d,*wherep,*noindfwpts,*wherefw,*match,*whoisfw,*nocoefs;     /* needed arrays */
  char *cltp;
  double *visc,*x[3];
  int *itrange; double *roundoff; /* use if available */
  int *coef_n, **coef_i; double **coef_c;  /* modify*/
  
  int its,ite,i,i4dp[4],iall,j,ip[4],iv[4],i4dm[4],iallm,ipt,id[4],ib,ic,ia[3],iptc[8],ibp,ibm,iw[8],nyok;
  int L,L2,L3, vtype=0,size;
  double tol=1.e-8;
  double vv[3][3],area[3],xmid[3],varea[3],xx[8][3],dx[3][3],grad[8][8][3],aa[3][3],vol[8];
  double dc,coefa[8][8];
  int two[3]={2,2,2};
  int dbug=0;
  double coeft[8][8],ppm,ppp;  /* for step (2-point) */
  /* for linear, geom8fx27 call and grad interp */
  double f[27][3],xf[27][3],fgrad[3],gradcen[8][3], fgc[8],fgcsum;  
  int iaf[3];
  
  nocoefs=(int *)need("nocoefs");
  
  if (jvisc>nocoefs[0]-1)
  { 
    printout("error coefvisc"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
    exitm4d(0);
  }
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallm=Prod4(i4dm);
  
  size=arraysize(namevisc);
  if (size==iallm) vtype=1;
  else if (size==6*iallm) vtype=6;
  else 
  { 
    printout("error coefvisc"," error, size of %s, %d, is not %d or %d (1 or 6x the between-the-pts array size)\n",namevisc,size,iallm,6*iallm); 
    exitm4d(0);
  }
  
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  visc=(double *)need(namevisc);
  
  x[0]=(double *)need("xyz"); 
  x[1]=x[0]+iall; 
  x[2]=x[1]+iall;
  coef_n=(int *)need("coef_n");
  coef_i=(int **)need("coef_i");
  coef_c=(double **)need("coef_c");
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  noindfwpts=(int *)need("noindfwpts");
  wherep=(int *)need("wherep");
  cltp=(char *)need("cltp");
  wherefw=(int *)need("wherefw");
  whoisfw=(int *)need("whoisfw");
  match=(int *)need("match");
  geom8init();  
  
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])  /* clear parts of array to be redone */
  if (coef_n[i]>0) {Loop(j,0,coef_n[i]) coef_c[i][j+coef_n[i]*jvisc]=0; }
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d) /* loop over all possible cont c.v. */
  {    
    iv[3]=ip[3];  id[3]=ip[3];
    if (dbug>0) printout("normal","ip %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]);
    if (wherep[In4(ip,i4dp)]<0) continue; /* skip if not valid cont c.v. */
    if (cltp[In4(ip,i4dp)]=='s' || cltp[In4(ip,i4dp)]=='S') continue;
    Loop(i,0,3) iv[i]=ip[i]-1;
    ipt=In4(iv,i4dm);
    Loop(i,0,3) Loop(j,0,3) vv[i][j]=0;
    if (vtype<2) Loop(i,0,3) vv[i][i]=visc[ipt]; /* isotropic */
    else   /* visc_i_j  given as 6 independent components 00 11 22 01 02 12 */ 
    {
      vv[0][0]=visc[ipt]; vv[1][1]=visc[ipt+iallm]; vv[2][2]=visc[ipt+2*iallm];
      vv[0][1]=visc[ipt+3*iallm]; vv[1][0]=vv[0][1];
      vv[0][2]=visc[ipt+4*iallm]; vv[2][0]=vv[0][2];
      vv[1][2]=visc[ipt+5*iallm]; vv[2][1]=vv[1][2];
    }
    if (dbug>0) printout("normal","vvn %lg %lg %lg\n",vv[0][0],vv[1][1],vv[2][2]);
    if (vv[0][0]+vv[1][1]+vv[2][2]==0) continue;    /* skip if no viscosity */
    Loop(i,0,8) Loop(j,0,8) coefa[i][j]=0;
    nyok=1;
    Loop3(ia,0,two)   /* grid at 8 corners */
    { 
      Loop(i,0,3) id[i]=iv[i]+ia[i];
      ib=ia[0]+2*ia[1]+4*ia[2];
      iptc[ib]=In4(id,i4d);
      iw[ib]=wherefw[iptc[ib]];
      if (iw[ib]<0) {nyok=0; break; }
      Loop(i,0,3) xx[ib][i]=x[i][iptc[ib]];
      if (dbug>1) printout("normal","ia %d %d %d ipt %d xyz %lg %lg %lg\n",
                         ia[0],ia[1],ia[2],iptc[ib],xx[ib][0],xx[ib][1],xx[ib][2]);
    }
    if (nyok==0) continue;  /* skip if some points have no equation */
    if (nylin==1) geom8fx27(f,xf,ip);  
    Loop(ib,0,8) iptc[ib]=match[iptc[ib]];
    Loop3(ia,0,two)   /* grad at 8 corners */
    {  
      ib=ia[0]+2*ia[1]+4*ia[2];
      Loop(i,0,8)  Loop(ic,0,3) grad[ib][i][ic]=0;
      Loop(i,0,3)   
      {
        dx[0][i]=xx[1+2*ia[1]+4*ia[2]][i]-xx[2*ia[1]+4*ia[2]][i];
        dx[1][i]=xx[ia[0]+2+4*ia[2]][i]-xx[ia[0]+4*ia[2]][i];
        dx[2][i]=xx[ia[0]+2*ia[1]+4][i]-xx[ia[0]+2*ia[1]][i];
      }
      Loop(i,0,3) {Cross(dx[(i+1)%3],dx[(i+2)%3],aa[i]); }
      vol[ib]=Dot(dx[0],aa[0]);
      if (dbug>1) printout("normal","ib %d vol %lg\n",ib,vol[ib]);
      if (dbug>1) {printout("normal","dx"); Loop(i,0,3) Loop(j,0,3) printout("normal"," %lg",dx[i][j]); printout("normal","\n"); }
      if (vol[ib]<=0) continue;  /* bad volume don't use assoc gradients */
      Loop(i,0,3) Loop(j,0,3) aa[i][j]/=vol[ib];
      Loop(i,0,3)
      { 
        grad[ib][1+2*ia[1]+4*ia[2]][i]+=aa[0][i];
        grad[ib][2*ia[1]+4*ia[2]][i]-=aa[0][i];
        grad[ib][ia[0]+2+4*ia[2]][i]+=aa[1][i];
        grad[ib][ia[0]+4*ia[2]][i]-=aa[1][i];
        grad[ib][ia[0]+2*ia[1]+4][i]+=aa[2][i];
        grad[ib][ia[0]+2*ia[1]][i]-=aa[2][i];
      }
    }
    if (nylin==0)
    {
      Loop(L,0,3)  /* each drection */
      { 
        L2=(L+1)%3; L3=(L+2)%3;
        id[L]=2*iv[L]+1;
        Loop(ia[L2],0,2) Loop(ia[L3],0,2)   /* each quarter area, ibp, ibm are points on either side*/
        { 
          ia[L]=1; ibp=ia[0]+2*ia[1]+4*ia[2];
          ia[L]=0; ibm=ia[0]+2*ia[1]+4*ia[2];
          if (dbug>1) printout("normal","L ia ib %d %d %d, ibp, ibm %d %d\n",L,ia[L2],ia[L3],ibp,ibm);
          if (match[iptc[ibp]]==match[iptc[ibm]]) continue;  /* skip if points use same eq */
          id[L2]=2*iv[L2]+ia[L2]; 
          id[L3]=2*iv[L3]+ia[L3];
          geom8areamid(area,xmid,id,L);
          Loop(i,0,3) varea[i]=Dot(vv[i],area);  /* viscosity times effective area */
          if (dbug>1) printout("normal","varea %lg %lg %lg\n",varea[0],varea[1],varea[2]);
          ia[L]=1; ibp=ia[0]+2*ia[1]+4*ia[2];
          ia[L]=0; ibm=ia[0]+2*ia[1]+4*ia[2];
          Loop(i,0,8)   /* coef using corner discretization separately for 2 points either side */
          { 
            coeft[ibm][i]=Dot(grad[ibm][i],varea);
            coeft[ibp][i]=Dot(grad[ibp][i],varea);
          }
          if (dbug>1) {printout("normal","coefm"); Loop(i,0,8) printout("normal"," %lg",coeft[ibm][i]); printout("normal","\n"); }
          if (dbug>1) {printout("normal","coefp"); Loop(i,0,8) printout("normal"," %lg",coeft[ibp][i]); printout("normal","\n"); }
          ppm=-coeft[ibm][ibm]+coeft[ibm][ibp]; 
          ppp=-coeft[ibp][ibm]+coeft[ibp][ibp];
          if (dbug>1) printout("normal"," ppm %lg ppp %lg\n",ppm,ppp);
          if (ppm<0) ppm=0;   if (ppp<0) ppp=0;
          if (ppm+ppp==0) continue;  /* omit as not being viscous */
          Loop(i,0,8)   /* favor the most stable corner discretization */
          { 
            dc=(ppm*coeft[ibm][i]+ppp*coeft[ibp][i])/(ppm+ppp);
            coefa[ibm][i]-=dc;
            coefa[ibp][i]+=dc;
          }
        }   /* each quarter area */
      } /* L */
    } /* nylin=0 */
    else /* nylin=1 */
    {
      Loop(L,0,3)  /* each drection */
      { 
        L2=(L+1)%3; L3=(L+2)%3;
        id[L]=2*iv[L]+1;
        Loop(ia[L2],0,2) Loop(ia[L3],0,2)   /* each quarter area, ibp, ibm are points on either side*/
        { 
          ia[L]=1; ibp=ia[0]+2*ia[1]+4*ia[2];
          ia[L]=0; ibm=ia[0]+2*ia[1]+4*ia[2];
          fgrad[L]=0; fgrad[L2]=0; fgrad[L3]=0;
          iaf[L]=1;
          for (iaf[L2]=ia[L2]; iaf[L2]<ia[L2]+2; iaf[L2]++)
            for (iaf[L3]=ia[L3]; iaf[L3]<ia[L3]+2; iaf[L3]++)
            { 
              ib=iaf[0]+3*iaf[1]+9*iaf[2];
              fgrad[L]+=.25*f[ib][L];
              fgrad[L2]+=.25*f[ib][L2];
              fgrad[L3]+=.25*f[ib][L3];
            }
          fgc[0]=(1-fgrad[0])*(1-fgrad[1])*(1-fgrad[2]);
          fgc[1]=fgrad[0]*(1-fgrad[1])*(1-fgrad[2]);
          fgc[2]=(1-fgrad[0])*fgrad[1]*(1-fgrad[2]);
          fgc[3]=fgrad[0]*fgrad[1]*(1-fgrad[2]);
          fgc[4]=(1-fgrad[0])*(1-fgrad[1])*fgrad[2];
          fgc[5]=fgrad[0]*(1-fgrad[1])*fgrad[2];
          fgc[6]=(1-fgrad[0])*fgrad[1]*fgrad[2];
          fgc[7]=fgrad[0]*fgrad[1]*fgrad[2];
          /* check for collapsed points */
          fgcsum=0;
          Loop(ic,0,8) 
          {
            if (vol[ic]<=0) fgc[ic]=0;
            fgcsum+=fgc[ic];
          }
          if (fgcsum==0) continue; /* shouldn't happen */
          Loop(ic,0,8) fgc[ic]/=fgcsum;
          Loop(ib,0,8) Loop(i,0,3)
          gradcen[ib][i]=fgc[0]*grad[0][ib][i]+fgc[1]*grad[1][ib][i]
          +fgc[2]*grad[2][ib][i]+fgc[3]*grad[3][ib][i]
          +fgc[4]*grad[4][ib][i]+fgc[5]*grad[5][ib][i]
          +fgc[6]*grad[6][ib][i]+fgc[7]*grad[7][ib][i];
          if (dbug>1) printout("normal","L ia ib %d %d %d, ibp, ibm %d %d\n",L,ia[L2],ia[L3],ibp,ibm);
          if (match[iptc[ibp]]==match[iptc[ibm]]) continue;  /* skip if points use same eq */
          id[L2]=2*iv[L2]+ia[L2]; 
          id[L3]=2*iv[L3]+ia[L3];
          geom8areamid(area,xmid,id,L);
          Loop(i,0,3) varea[i]=Dot(vv[i],area);  /* viscosity times effective area */
          if (dbug>1) printout("normal","varea %lg %lg %lg\n",varea[0],varea[1],varea[2]);
          ia[L]=1; ibp=ia[0]+2*ia[1]+4*ia[2];
          ia[L]=0; ibm=ia[0]+2*ia[1]+4*ia[2];
          Loop(i,0,8)    /* using grad linearly interpolated to area center */
          {
            dc=Dot(gradcen[i],varea);
            coefa[ibm][i]-=dc;
            coefa[ibp][i]+=dc;
          }
        }   /* each quarter area */
      } /* L */
    } /* nylin=1 */
    if (dbug>0) 
    {
      printout("normal","coefa for ipt"); Loop(i,0,8) printout("normal"," %d",iptc[i]); printout("normal","\n");
      Loop(i,0,8) {Loop(j,0,8) printout("normal"," %lg",coefa[i][j]); printout("normal","\n"); }
    }
    Loop(i,0,8) coefcadd(iw[i],coef_n,coef_i,coef_c,nocoefs[0],jvisc,8,iptc,coefa[i]);
  } /* end 	ip loop over cont. c.v. */	  
  
  Loop(i,noindfwpts[its+1],noindfwpts[ite+2])    /* clean up redundent coefficients */
  {
    coefcombine(i,coef_n,coef_i,coef_c,nocoefs[0],tol);
    coeforder(whoisfw[i],i,coef_n,coef_i,coef_c,nocoefs[0]);  
  }
}

