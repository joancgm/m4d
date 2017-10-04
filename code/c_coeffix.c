/* contains c_coeffix */
#include "global.h"
/* calc sum pose and disp of coefs at jfrom, and set fixed coefs at  jto (if jto>=0)
 Note:  set jto=0 if eqnsolve is to solve with these coefficients  */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

void c_coeffix(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*match,*noindfwpts,*whoisfw,*nocoefs,*wherefw;  /* needed arrays */
  char *clt;
  int *itrange; double *roundoff; /* use if available */
  int *cn, **ci;  double **cc;  /* need and modify  if jto>=0 */
  double *sum,*pose,*disp;  /* calculate */
  int jfrom,jto,limsumnz;  /* input parameters */
  double poselim,displim;
  int nomit; char omit[20];
  double cadd[2]; int iptc[2];
  
  int its,ite;
  int i,j,k,iall,ipt,iw,itot,ipose[11],itotsum,itotpsum,itotpose,itotdisp,iwa;
  double rr=1.e-10;
  double cen,cp,sump,cmax,cmin,sumamax,posemin,dispmin,dcen,dcend,dcenp,ch;
  int bug=0; int nyperm=1;  /* no longer an input */
  double fc,dcensum;
  
  jfrom=readint(fpin);
  jto=readint(fpin);
  limsumnz=readint(fpin);
  poselim=readdouble(fpin);
  displim=readdouble(fpin);
  poselim=min(poselim,1);
  displim=min(displim,-1);
  
  printout("normal"," nyperm %d, coef %d to %d. ",nyperm,jfrom,jto);
  if (jto==jfrom) printout("normal"," warning with jto=jfrom, mixing might be direction biased\n");
  printout("normal"," Limits: limsumnz %d  poselim %lg displim %lg\n",limsumnz,poselim,displim);
  nomit=readint(fpin);
  printout("normal","omit %d pt types:",nomit);
  Loop(i,0,nomit) {omit[i]=read1charname(fpin); printout("normal"," %c",omit[i]); }
  
  nocoefs=(int *)need("nocoefs");
  if (max(jfrom,jto)>nocoefs[0]-1)
  { 
    printout("error c_coeffix"," error, coef only dimensioned for %d coefs, change with coefinit\n",nocoefs[0]);
    exitm4d(0);
  }
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  itrange=(int *)find("itrange");
  if (itrange==0) {its=0; ite=i4d[3]-1; }
  else {its=itrange[0]; ite=itrange[1]; }
  printout("normal"," ,  its %d ite %d\n", its,ite);
  
  noindfwpts=(int *)need("noindfwpts");
  match=(int *)need("match");
  whoisfw=(int *)need("whoisfw");
  wherefw=(int *)need("wherefw");
  cn=(int *)need("coef_n");
  ci=(int **)need("coef_i");
  cc=(double **)need("coef_c");
  roundoff=(double *)find("roundoff");
  if (roundoff>0) rr=roundoff[0];
  clt=(char *)need("clt");
  
  if (nyperm>0)
  {
    sum=(double *) createarray("sum",iall,'d',0);
    pose=(double *) createarray("pose",iall,'d',0);
    disp=(double *) createarray("disp",iall,'d',0);
  }
  else 
  {
    sum=(double *)tmalloca(iall,'d');
    pose=(double *)tmalloca(iall,'d');
    disp=(double *)tmalloca(iall,'d');
  }
  Loop(i,0,iall) {sum[i]=0; pose[i]=1; disp[i]=0;}  /* init no equation */
  
  /* create center point coefficient if not there and make sure cener is first one 
   then orcer coefs by value  and set jto ones =0  if appropriate */
  Loop(iw,noindfwpts[its+1],noindfwpts[ite+2])   
  { 
    if (cn[iw]>0)
    { 
      ipt=whoisfw[iw]; 
		coefcenterfirst(ipt,iw,cn,ci,cc,nocoefs[0]);
		coeforderv(ipt,iw,cn,ci,cc,nocoefs[0],jfrom);
		if (jto>=0 && jfrom!=jto)
        Loop(k,0,cn[iw]) cc[iw][k+jto*cn[iw]]=0;
    }
  }
  
  itotsum=0; itotpsum=0; sumamax=0; 
  itot=0; posemin=2; dispmin=0;
  itotpose=0; itotdisp=0;
  Loop(k,0,11) ipose[k]=0;
  /*   calc sums and fix first before doing pose,disp, all points*/
  Loop(iw,noindfwpts[its+1],noindfwpts[ite+2])  
  { 
    if (cn[iw]==0) continue; 
    ipt=whoisfw[iw];
    /* calc sum and fix if appropriate */
    sump=0;
    Loop(k,0,cn[iw])
    { 
      cp=cc[iw][k+cn[iw]*jfrom];
      sum[ipt]+=cp;
      if (cp>0) sump+=cp;
    }
    if (abs(sum[ipt])<rr*sump) sum[ipt]=0;
    dcensum=0;
    if (sum[ipt]!=0) 
    {
      if ( jto>=0 && (sum[ipt]<0 || limsumnz>0)) 
		{  
        cc[iw][jto*cn[iw]]-=sum[ipt];
        dcensum=-sum[ipt];
        if (sum[ipt]>0) itotpsum++;
        else itotsum++;
		}
      sum[ipt]/=(sump+abs(sum[ipt]));
      if (abs(sum[ipt])> sumamax) sumamax=abs(sum[ipt]);
    }
    /* pose and disp calcs */
    if (cn[iw]<2) continue; 
    j=0;  /* don't check for specified point types */
    Loop(k,0,nomit) if (clt[ipt]==omit[k]) j=1;
    if (j==1) continue;
    sump=0; cmax=0; cmin=0;
    cen=cc[iw][cn[iw]*jfrom];  if (jto!=jfrom) cen +=dcensum;
    if (cen>0) { sump=cen; cmax=cen;}
    if (cc[iw][1+cn[iw]*jfrom]>cmax) cmax=cc[iw][1+cn[iw]*jfrom];
    cmin=cc[iw][cn[iw]-1+cn[iw]*jfrom];
    if (cmax==cmin) continue;
    Loop(k,1,cn[iw]) 
    { 
      cp=cc[iw][k+cn[iw]*jfrom];
      if (cp<=0) break;
      sump+=cp;
    }
    if (sump<=0) {pose[ipt]=-10; disp[ipt]=-10; }
    else
    { 
      pose[ipt]=cen/sump;
      if (cen>0) disp[ipt]=cmin/cen;
      else disp[ipt]=-10;
    }
    if (jto<0 && pose[ipt]<0) printout("normal","warning: ipt %d pose %lg disp %lg\n",ipt,pose[ipt],disp[ipt]);
    itot++;
    posemin=min(pose[ipt],posemin);
    dispmin=min(disp[ipt],dispmin);
    Loop(k,0,11) if (pose[ipt]<(double)(0.1*k)) {ipose[k]++; break;}
    if (pose[ipt]>=poselim && disp[ipt]>=displim) continue;
    /* pose and disp calcs */
    if (bug>0) 
    { 
      int iexp[4];
      iexpand(ipt,i4d,iexp);
      printout("normal","ipt %d %d %d %d pose %lg disp %lg  cen %lg sum %lg\n",
              iexp[0],iexp[1],iexp[2],iexp[3],pose[ipt],disp[ipt],cen,sum[ipt]);
      Loop(k,0,cn[iw]) 
      {
        iexpand(ci[iw][k],i4d,iexp);
        if (cc[iw][k+cn[iw]*jfrom]>0 || cc[iw][k+cn[iw]*jfrom]< -cen) 
          printout("normal"," ****");
        else  printout("normal", "          ");
        printout("normal"," %d %d %d %d %lg %lg\n",
                iexp[0],iexp[1],iexp[2],iexp[3],cc[iw][k+cn[iw]*jfrom],cc[iw][k+cn[iw]*jto]);
      }
    } /* bug */
    if (jto>=0)  /* calculate dcen for needed change to center point */
    { 
      if (bug>0) printout("normal","cmin %lg cen %lg sump %lg displim %lg poselim %lg\n",
                         cmin,cen,sump,displim,poselim);
      dcend=cmin/displim-cen;
      dcenp=sump*poselim-cen;
      dcen=max(dcend,dcenp);
      if (bug>0) printout("normal","dcend %lg dcenp %lg dcen %lg\n",dcend,dcenp,dcen); 
      if (dcen<=0) continue;
      if (dcenp>dcend) itotpose++;  else itotdisp++;
      fc=min(dcen/(sump-max(cen,0)),1);
      if (bug>0) printout("normal"," fc %lg\n",fc);
      Loop(k,1,cn[iw])
      { 
        cp=cc[iw][k+cn[iw]*jfrom];
        if (cp<rr*sump) break;
        ch=fc*cp;
        cadd[0]=ch;
        cadd[1]=-ch;
        iptc[0]=ci[iw][0];
        iptc[1]=ci[iw][k];
        iwa=wherefw[ci[iw][k]];
        coefcadd(iw,cn,ci,cc,nocoefs[0],jto,2,iptc,cadd);
        cadd[0]=-ch;
        cadd[1]=ch;
        coefcadd(iwa,cn,ci,cc,nocoefs[0],jto,2,iptc,cadd);
      }
    } /* done dcc for fixed  coefs */
  } /* iw loop */
  
  if (sumamax>0) printout("normal","abs(sum)max %lg, fixed+ at %d points, fixed- at %d points\n",
                         sumamax,itotpsum,itotsum);
  printout("normal","for %d points, pose min %lg disp min %lg\n", itot,posemin,dispmin);
  printout("normal","pose groups");
  Loop(k,0,11) printout("normal","  %d <%lg,",ipose[k],0.1*k);
  printout("normal","\n");
  if (itotpose>0 || itotdisp>0)
    printout("normal","fixed %d eqs for pose, %d eqs for disp\n",
            itotpose,itotdisp);
  if (nyperm>0)
    Loop(i,0,iall)  /* make arrays print and plottable */
  {
    sum[i]=sum[match[i]]; 
    pose[i]=pose[match[i]]; 
    disp[i]=disp[match[i]]; 
  }
  if (nyperm<=0) tmalloca(-1,'d'); 
  
  if (jto>=0 && jto!=jfrom)  /* add jfrom to jto */
    Loop(iw,noindfwpts[its+1],noindfwpts[ite+2])   
  { 
    if (cn[iw]>0)
		Loop(k,0,cn[iw]) cc[iw][k+jto*cn[iw]]+=cc[iw][k+jfrom*cn[iw]];
  }
}

