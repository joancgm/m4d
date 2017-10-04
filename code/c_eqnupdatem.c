/* contains c_eqnupdatem */
#include "global.h"
#include "eqnsubs.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/* a single center point update of the velocities 
 equivalent to eqnsolves with procedure c 1 and then adding the changes 
 to the velocity 
 csym used if available for symmetry planes
 read omit to add point types to be omitted
 any point which lack coefficients will also be left unchanged 
 */

void c_eqnupdatem(FILE *fpin, FILE *fprint)
{  
  int *i4d,*npts,*neqs,*where,*whois,*cn,**ci, *match; /* needed arrays */
  double **cc, *u, *rhsin; 
  char *clt,*csym;  
  double *ch=0; /* create or temporary */
  int *itrange; double *camin; /* use if  available */
  
  int  nysavch,nomit;  /* input parameters */
  char omit[20];
  
  double *rhs,*cplus,*csave=0; /* temporary arrays */
  
  int i,j,k,L,ieqstart,ieqend,its,ite,itot,itota;
  int ii[4],is[4],ie[4];
  char *ny=0;
  int imin,imax,iimin[4],iimax[4];
  double err[3],errinit[3];
  char name[3]="U0",namerhs[6]="rhsU0",namedu[4]="dU0",c123[4]="123",cxyz[4]="xyz";
  
  nysavch=readint(fpin);     /* read input parameters */
  nomit=readint(fpin);
  printout("normal","ny save dUn %d omit %d pt types:",nysavch,nomit);
  Loop(i,0,nomit) {omit[i]=read1charname(fpin); printout("normal"," %c",omit[i]); }
  printout("normal","\n");
  
  npts=(int *)need("nogpts");   /* get needed arrays */
  i4d=(int *)need("idim4d");
  neqs=(int*)need("noindfwpts");
  where=(int *)need("wherefw");
  whois=(int *)need("whoisfw");
  match=(int *)need("match");
  clt=(char *)need("clt");
  csym=(char *)need("csym");
  cn=(int *)need("coef_n");
  ci=(int **)need("coef_i");
  cc=(double **)need("coef_c");
  camin=(double *)find("cam");
  itrange=(int *)find("itrange");
  its=0; ite=i4d[3]-1;
  if (itrange>0) {its=itrange[0]; ite=itrange[1]; }
  ieqstart=neqs[its+1]; 
  ieqend=neqs[ite+2];
  
  if (camin>0)   /* add camin to center coef */
  {  
    csave=(double*)tmalloca(neqs[0],'d');
    Loop(i,ieqstart,ieqend)
    if (cn[i]>0)
    { 
      csave[i]=cc[i][0];
      cc[i][0]+=camin[i];
    }
  }
  else printout("normal","no cam center pt additions\n");
  cplus=(double *)tmalloca(neqs[0],'d');
  eqncplus(ieqstart,ieqend,cn,cc, cplus); /* cplus=sum of + coef */
  rhs=(double *)tmalloca(neqs[0],'d');
  if (nysavch<=0) ch=(double *)tmalloca(npts[0],'d');
  Loop(L,0,3)  /* for U1, U2, U3 momentum equations */
  {  
    name[1]=c123[L]; namerhs[4]=c123[L]; namedu[2]=c123[L];
    
    if (i4d[L]==2 && csym[2*L]==cxyz[L] && csym[2*L+1]==cxyz[L])
    { 
      printout("normal","2-d symmetry plane, omit update for %s\n",name);
      continue;
    }
    u=(double *)need(name);
    rhsin=(double *)need(namerhs);
    Loop(i,ieqstart,ieqend) rhs[i]=rhsin[i];
    if (nysavch>0) ch=(double *)createarray(namedu,npts[0],'d',0);
    Loop(i,0,npts[0]) ch[i]=0;
    
    itot=0; itota=0;
    if (nomit>0) /* change sign of cplus so that eqs are omitted */
    {  
      Loop(i,ieqstart,ieqend)
      {  
        itota++;
        if (cplus[i]>0)
          Loop(j,0,nomit) 
          if (clt[ci[i][0]]==omit[j]) 
          {cplus[i]=-cplus[i]; itot++; break;}
      }	
    }
    /* check symmetry */
    Loop(i,0,3) Loop(j,0,2)
    if (csym[j+2*i]==cxyz[L])
    { 
      Loop(k,0,3) {is[k]=0; ie[k]=i4d[k];} is[3]=its; ie[3]=ite+1;
      if (j==0) ie[i]=1;
      else is[i]=ie[i]-1;
      Loop(ii[3],is[3],ie[3]) Loop(ii[2],is[2],ie[2])
      Loop(ii[1],is[1],ie[1]) Loop(ii[0],is[0],ie[0])
      { 
        k=where[In4(ii,i4d)];
        if (k>=0) if (cplus[k]>0) {cplus[k]=-cplus[k]; itot++; }
      }
    }
    
    printout("normal"," update %s, %d of %d eqs, type and sym-plane deactivated\n",name,itot,itota);
    
    eqnerrorwho(ieqstart,ieqend,rhs,cplus,ny, err,&imin,&imax); /* error analysis */ 
    iexpand(whois[imin],i4d,iimin); iexpand(whois[imax],i4d,iimax);
    printout("normal"," init:  err-rms %lg    err-min %lg at %d %d %d %d    err-max %lg at %d %d %d %d\n",
            err[0],err[1],iimin[0],iimin[1],iimin[2],iimin[3],
            err[2],iimax[0],iimax[1],iimax[2],iimax[3]);
    Loop(i,0,3) errinit[i]=err[i];
    
    eqncenter(ieqstart,ieqend,cplus,rhs,ci, ch);
    eqnrhs(ieqstart,ieqend,cn,ci,cc,ch,rhsin,  rhs); 
    
    eqnerrorwho(ieqstart,ieqend,rhs,cplus,ny, err,&imin,&imax); /* error analysis */ 
    iexpand(whois[imin],i4d,iimin); iexpand(whois[imax],i4d,iimax);
    printout("normal"," center update  err-rms %lg    err-min %lg at %d %d %d %d    err-max %lg at %d %d %d %d\n",
            err[0],err[1],iimin[0],iimin[1],iimin[2],iimin[3],
            err[2],iimax[0],iimax[1],iimax[2],iimax[3]);
    
    Loop(i,ieqstart,ieqend) if (cn[i]>0) u[ci[i][0]]+=ch[ci[i][0]]; /* update */
    Loop(i,0,npts[0]) u[i]=u[match[i]]; /* set matched points */
    if (itot>0 && L<2)
      Loop(i,ieqstart,ieqend) if (cplus[i]<0) cplus[i]=-cplus[i];
  }
  
  /* end of solve and update , clean up*/
  if (camin>0)   /* reinstate center coef */
    Loop(i,ieqstart,ieqend)
    if (cn[i]>0) cc[i][0]=csave[i];
}

