/* contains c_contcpcdu */
#include "global.h"

#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
/* continuity equations, pressure-correction coefficients  for between the points volumes */
/* initial setting of cpc_n,i,c   based on du : c0 on pts du; c1 full midface delp  correction */
/* calc  coefs for pc based on abrev momentum */

void c_contcpcdu(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*wherep,*noindppts,*matchpc; char *csym, *clt;	 	/* needed arrays */ 
  double *cam,*cplus,*rho,*x[3],*xyzp,*volmom;
  int *wherefw,*cpflop_n, **cpflop_i; double **cpflop_c;
  int *itrange; double *roundoff; /* use if available */
  int *cpc_n, **cpc_i; double **cpc_c;    /* create and/ or update */
  
  int i4dp[4],ip[4],i,j,iall,ieq,L,L2,L3,ig[4],ia[4],sign,ib,ic,inew,n,nfpt;
  int mpc,ipa[4],mpca,iallp,nst[8],nadd,kpc=-1,kpca=-1,nyrc,iptc,iptca;
  double area[4][3],da[3],dxpc[3],dxpcsq=0,flopc[3],radcam[3],facrc=0,cenadd;
  double tol=1.e-8;
  int its,ite;
  int db=0;   /* reset to 1 to get debug prints */
  
  i4d=(int *)need("idim4d");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iall=Prod4(i4d);
  iallp=Prod4(i4dp);
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  rho=(double *)need("rho");
  geomcinit();  /* for areas */
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall; x[2]=x[0]+2*iall;
  csym=(char *)need("csym");
  clt=(char *)need("clt");
  wherep=(int *)need("wherep");
  noindppts=(int *)need("noindppts");
  
  matchpc=(int *)need("matchpc");
  xyzp=(double *)need("xyzp");
  volmom=(double *)need("volmom");
  cam=(double *)need("cam");
  cplus=(double *)need("cplus");
  wherefw=(int *)need("wherefw"); 
  cpflop_n=(int*)need("cpflop_n"); 
  cpflop_i=(int**)need("cpflop_i"); 
  cpflop_c=(double**)need("cpflop_c"); 
  
  cpc_n=(int *)find("cpc_n");      
  n=noindppts[0];
  if (cpc_n==0)    /* create if doesn't exist */
  { 
    cpc_n=(int *)createarray("cpc_n",n,'i',0);
    cpc_i=(int **)createarray("cpc_i",n,'p',0);
    cpc_c=(double **)createarray("cpc_c",n,'p',0);
    Loop(i,0,n) {cpc_n[i]=0; cpc_i[i]=0; cpc_c[i]=0;}
  }
  else
  { 
    cpc_i=(int**)need("cpc_i");   
    cpc_c=(double**)need("cpc_c"); 
    Loop(ieq,noindppts[its+1],noindppts[ite+2])
    if (cpc_n[ieq]>0) 
    {  
      free(cpc_i[ieq]); free(cpc_c[ieq]); 
      cpc_n[ieq]=0; cpc_i[ieq]=0; cpc_c[ieq]=0;
    }
  }
  
  Loop(ip[3],its,ite+1) Loop3(ip,1,i4d)    /* -----------over all cont c.v. --------------*/
  { 
    mpc=In4(ip,i4dp);
    ieq=wherep[mpc];
    if (db>0)   printout("normal"," ip %d %d %d ieq %d mpc %d cpcn %d\n",
                       ip[0],ip[1],ip[2],ieq,mpc,cpc_n[ieq]);
    if (ieq>=0) /* do only for active volumes */
    { 
      inew=0; 
      ig[3]=ip[3];      /* determine enough space for the coefficients */
      Loop(ia[2],0,2) Loop(ia[1],0,2) Loop(ia[0],0,2)
      {
        Loop(i,0,3) ig[i]=ip[i]-1+ia[i];   /* momentum point */
        n=wherefw[In4(ig,i4d)];
        nst[ia[0]+2*ia[1]+4*ia[2]]=inew;
        if (n>=0) inew+=cpflop_n[n];
      }
      cpc_n[ieq]=inew;
      if (inew>0)  /* allocate space and initialize */
      { 
        cpc_i[ieq]=(int *)smalloca(inew,'i');
        cpc_c[ieq]=(double *)smalloca(inew*2,'d');
        Loop(i,0,2*inew) cpc_c[ieq][i]=0;
        Loop(ia[2],0,2) Loop(ia[1],0,2) Loop(ia[0],0,2)
        {
          Loop(i,0,3) ig[i]=ip[i]-1+ia[i];   
          n=wherefw[In4(ig,i4d)];
          if (n>=0) Loop(i,0,cpflop_n[n]) 
            cpc_i[ieq][nst[ia[0]+2*ia[1]+4*ia[2]]+i]=cpflop_i[n][i];
        }
        /* }  defer ending inew to the end */
        if (db>0) printout("normal","new %d nst1 nst8 %d %d %d %d %d %d %d %d\n",
                         inew,nst[0],nst[1],nst[2],nst[3],nst[4],nst[5],nst[6],nst[7]);
        
        Loop(L,0,3)   /* each direction */
        { 
          L2=(L+1)%3, L3=(L+2)%3;
          Loop(ia[L],0,2) /* each surface */
          {  
            sign=2*ia[L]-1;  
            ig[L]=ip[L]-1+ia[L];
            /* ------  check if can use mid-face (rc) corrections ------*/
            Loop(i,0,4) ipa[i]=ip[i];    /* find pc location on other side */
            ipa[L]+=sign;
            mpca=In4(ipa,i4dp);
            if (db>0) printout("normal","L %d ia[L] %d mpca %d new %d \n",L,ia[L],mpca,inew);
            nyrc=1;   /* yes to start for ny Rhie_Chow coefs as well */
            if (matchpc[mpca]>=0) /* if sleeping point, find real point on other side */
            { 
              if (matchpc[mpca]==mpc) mpca=matchpc[mpca+iallp];
              else if (matchpc[mpca+iallp]==mpc) mpca=matchpc[mpca];
              else nyrc=0;
            }
            if (nyrc==1)  /* find cpc indices to be  used for the rc corrections  */
            {
              kpc=-1; kpca=-1;
              Loop(i,0,inew) { if (cpc_i[ieq][i]==mpc) {kpc=i; break;}}
              Loop(i,0,inew) {if (cpc_i[ieq][i]==mpca) {kpca=i; break;}}
              if (kpc==-1 || kpca==-1 || kpc==kpca) nyrc=0;
            }
            if (nyrc==1)    
            { 
              Loop(i,0,3) dxpc[i]=sign*(xyzp[mpca+i*iallp]-xyzp[mpc+i*iallp]);
              if ((ig[L]==0 || ig[L]==i4d[L]-1)
                  && (csym[6+2*L]=='r' || csym[6+2*L]=='R'))
              { 
                ig[L2]=ip[L2]-1; ig[L3]=ip[L3]-1; iptc=In4(ig,i4d);
                ig[L]=i4d[L]-1-ig[L]; iptca=In4(ig,i4d);
                Loop(i,0,3) dxpc[i]=sign*
                (xyzp[mpca+i*iallp]-x[i][iptca]-xyzp[mpc+i*iallp]+x[i][iptc]);
                if (db>0) printout("normal","repeat dxpc %lg %lg %lg\n",dxpc[0],dxpc[1],dxpc[2]);
              }
              dxpcsq=dxpc[0]*dxpc[0]+dxpc[1]*dxpc[1]+dxpc[2]*dxpc[2];
              if (dxpcsq==0) nyrc=0;
            }
            if (db>0) printout("normal","nyrc %d mpca %d kpc %d kpca %d\n",nyrc,mpca,kpc,kpca);
            /*-------- end rc checks-------*/
            Loop(i,0,3) ig[i]=ip[i]-1;     /* get areas */
            ig[L]=ig[L]+ia[L];
            geomcavector(area,ig,L);
            
            Loop(ia[L2],0,2) Loop(ia[L3],0,2) /* indices of 4 corners of surface */
            { 
              Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
              iptc=In4(ig,i4d);
              ib=ia[L2]; ic=ia[L3];
              nfpt=wherefw[iptc];  
              if (nfpt>=0 && clt[iptc] !='w')  /* make sure there is a eq for it, omit wall pts */
              {  
                Loop(i,0,3) da[i]=-area[ib+2*ic][i]*sign; /* now in/out signed */
                if (cplus[nfpt]>0) /* for rho da/cam */
                  Loop(i,0,3) radcam[i]=rho[iptc]*da[i]/(cam[nfpt]+cplus[nfpt]);
                else Loop(i,0,3)  radcam[i]=0;
                nadd=nst[ia[0]+2*ia[1]+4*ia[2]];
                n=cpflop_n[nfpt];
                
                if (nyrc==1)   /* part 2 rc coefficients */
                { 
                  facrc=Dot(dxpc,radcam)/dxpcsq;
                  cenadd=-facrc*volmom[nfpt]*sign;
                  if (db>0) printout("normal","facrc vol cenadd %lg, %lg %lg\n",
                                   facrc,volmom[nfpt],cenadd);
                  cpc_c[ieq][kpc+inew]+=cenadd;
                  cpc_c[ieq][kpca+inew]-=cenadd;
                }
                
                Loop(j,0,n)   /* each cpflop coef */
                { 
                  Loop(i,0,3) flopc[i]=cpflop_c[nfpt][j+n*i];
                  cpc_c[ieq][j+nadd]+=Dot(flopc,radcam);
                  if (nyrc==1) /* set part 2 rc coefficients */
                    cpc_c[ieq][j+nadd+inew]-=facrc*Dot(dxpc,flopc);
                } /* each cpflop coef */
              }  /* nfpt check */
            } /* each corner point */
            
          } /* each surface */
        } /* each direction */
        
        if (db>0) 
        { 
          printout("normal","for p pt %d %d %d, ieq %d no of coef %d:   i j k n c\n",
                 ip[0],ip[1],ip[2],ieq,cpc_n[ieq]);
          Loop(i,0,cpc_n[ieq])
          { 
            iexpand(cpc_i[ieq][i],i4dp,ia);
            printout("normal","        %d    %d %d %d %d     %lg  %lg    sum   %lg\n",
                   cpc_i[ieq][i],ia[0],ia[1],ia[2],ia[3],cpc_c[ieq][i],cpc_c[ieq][i+inew],
                   cpc_c[ieq][i]+cpc_c[ieq][i+inew]);
          }
        }
        coefcombine(ieq,cpc_n,cpc_i,cpc_c,2,tol);
        coeforder(mpc,ieq,cpc_n,cpc_i,cpc_c,2); 
        if (db>0) 
        { 
          printout("normal","for p pt %d %d %d, ieq %d no of coef %d:   i j k n c\n",
                 ip[0],ip[1],ip[2],ieq,cpc_n[ieq]);
          Loop(i,0,cpc_n[ieq])
          { 
            iexpand(cpc_i[ieq][i],i4dp,ia);
            printout("normal","        %d    %d %d %d %d     %lg  %lg  sum   %lg\n",
                   cpc_i[ieq][i],ia[0],ia[1],ia[2],ia[3], cpc_c[ieq][i],  cpc_c[ieq][i+inew], 
                   cpc_c[ieq][i]+cpc_c[ieq][i+inew]);
          }
        }
      }   /* inew>0) */
    } /* active volume */
  }  /* ip loop */
}
