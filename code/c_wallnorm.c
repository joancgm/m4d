/* contains x3line, c_wallnorm */
/* determine wall point list and other items needed to set prop bndry conditions */
#include "global.h"

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b)
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

/*----x3line--xx=interpolated paralled direction at middle of 3 pts--------*/
void x3line(int m1, int m2, int m3, int iall, double *xyz, double *xx)
{ 
  int n,k;
  double x[2][3],r[2];
  for (n=0;n<3;n++)
  {
    x[0][n]=xyz[m2+n*iall]-xyz[m1+n*iall];
    x[1][n]=xyz[m3+n*iall]-xyz[m2+n*iall];
  }
  for (k=0;k<2;k++) r[k]=sqrt(x[k][0]*x[k][0]+x[k][1]*x[k][1]+x[k][2]*x[k][2]);
  if (r[0]==0) for (n=0;n<3;n++) xx[n]=x[1][n];
  else if (r[1]==0) for (n=0;n<3;n++) xx[n]=x[0][n];
  else for (n=0;n<3;n++) xx[n]=x[0][n]*r[1]/r[0]+x[1][n]*r[0]/r[1];
}

/*----------------------------------------------*/
void c_wallnorm(FILE *fpin, FILE *fprint)
{ 
  int *i4d,*match; char *csym; char *clt; double *xyz;  /* need */
  int *noindwpts, *whoisw, *wnear; double *wnorm;   /* set */
  
  int iall3,iall,i,j=0,k,m,n,L,ifac[3],ii[4],inear[4];
  int mm[3][2],nyw[3],nywt,nyf[3],nyft,notdone=0,nychecksym=0; 
  int ia,ib,ic[3],L2,L3,idif[27],iflow[27];
  char cnear[3][2],cxyz[4]="xyz"; 
  double dx[3][3],wmag;
  int bug=0;
  
  i4d=(int *)need("idim4d");
  iall3=i4d[0]*i4d[1]*i4d[2];
  iall=iall3*i4d[3];
  ifac[0]=1; ifac[1]=i4d[0]; ifac[2]=i4d[0]*i4d[1];
  match=(int *)need("match");
  clt=(char *)need("clt");
  xyz=(double *)need("xyz");
  csym=(char *)need("csym");
  Loop(i,0,6) Loop(m,0,3) if (csym[i]==cxyz[m]) nychecksym=1;
  
  /* noindwpts and whoisw */
  noindwpts=(int *)createarray("noindwpts",i4d[3]+2,'i',0);
  Loop(i,0,i4d[3]+2) noindwpts[i]=0;
  Loop(i,0,iall) if (match[i]==i && clt[i]=='w') 
  {noindwpts[0]++; noindwpts[2+i/iall3]++;}
  printout("normal","nunber of independent wall point = %d\n",noindwpts[0]);
  
  whoisw=(int *)createarray("whoisw",noindwpts[0],'i',0);
  k=0;
  Loop(i,0,iall)  if (match[i]==i && clt[i]=='w')
  {whoisw[k]=i; k++;}
  
  /* wall normal vectors */
  wnorm=(double *)createarray("wnorm",3*noindwpts[0],'d',0);
  wnear=(int *)createarray("wnear",noindwpts[0],'i',0);
  
  Loop(i,0,noindwpts[0])
  {
    m=whoisw[i];
    wnorm[3*i]=0; wnorm[3*i+1]=0; wnorm[3*i+2]=0; wnear[i]=-1;
    iexpand(m,i4d,ii);
    if (bug>1) printout("normal","wall pt %d %d %d %d\n",ii[0],ii[1],ii[2],ii[3]);
    /* find  neighbor indices each direction mm[3][2], -1 if outside
     and point type, cnear[3][2] s,w,f,o(outside) */
    nywt=0; nyft=0;
    Loop(L,0,3)
    {
      nyw[L]=0; dx[L][0]=0; dx[L][1]=0; dx[L][2]=0;
      /* set - side */
      if (ii[L]==0) {mm[L][0]=-1; cnear[L][0]='o'; }
      else 
      {
        mm[L][0]=match[m-ifac[L]]; 
        if (clt[mm[L][0]]=='s') cnear[L][0]='s';
        else if (clt[mm[L][0]]=='w') cnear[L][0]='w';
        else cnear[L][0]='f';
      }
      /* set + side */
      if(ii[L]==i4d[L]-1) {mm[L][1]=-1; cnear[L][1]='o'; }
      else 
      { 
        if (csym[6+2*L]=='R') /* leave in 2-d to get direction */
        {
          mm[L][1]=m+ifac[L]; cnear[L][1]='w'; 
        }
        else 
        {
          mm[L][1]=m; j=0;
          while(mm[L][1]>=0 && match[mm[L][1]]==m)
          { 
            j++; 
            if (ii[L]+j==i4d[L]-1)  {mm[L][1]=-1; cnear[L][1]='o'; }
            else mm[L][1]+=ifac[L];
          }
          if (mm[L][1]>0)
          {
            mm[L][1]=match[mm[L][1]];
            if (clt[mm[L][1]]=='s') cnear[L][1]='s';
            else if (clt[mm[L][1]]=='w') cnear[L][1]='w';
            else cnear[L][1]='f';
          }
        }
      }
      nyw[L]=1; nywt++;    /* create dx line segments */
      if (cnear[L][0]=='w' && cnear[L][1]=='w')
        x3line(mm[L][0],m,mm[L][1],iall,xyz,dx[L]);
      else if (cnear[L][0]=='w')
        Loop(j,0,3) dx[L][j]=xyz[m+iall*j]-xyz[mm[L][0]+iall*j];
      else if (cnear[L][1]=='w')
        Loop(j,0,3) dx[L][j]=xyz[mm[L][1]+iall*j]-xyz[m+iall*j];
      else {nyw[L]=0; nywt--;}
      nyf[L]=-1;
      if (cnear[L][0]=='f') {nyf[L]=0; nyft++; }
      else if (cnear[L][1]=='f') {nyf[L]=1; nyft++; }
      if (bug>1) printout("normal","L %d mm clt %d %c %d %c dx %lg %lg %lg\n",
                         L,mm[L][0],cnear[L][0],mm[L][1],cnear[L][1],dx[L][0],dx[L][1],dx[L][2]);
    } /* L */
    if (nywt==3) 
    { 
      if (nyft==3) /* corner, find quadrant and averave vectors */
		{ 
        if (bug>0) printout("normal","warning nywt %d nyft %d not done\n",nywt,nyft);
		}
		else if (nyft==2) /* edge form new 3-pt line and cross */
		{ 
        Loop(L,0,3) if (cnear[L][0]!='f' && cnear[L][1]!='f') break;
        L2=(L+1)%3; L3=(L+2)%3;
        if (cnear[L2][0]=='f') ia=0; else ia=1;
        if (cnear[L3][0]=='f') ib=0; else ib=1;
        wnear[i]=match[m+(2*ia-1)*ifac[L2]+(2*ib-1)*ifac[L3]];
        if (ia==ib) x3line(mm[L2][1-ia],m,mm[L3][1-ib],iall,xyz,dx[L2]);
        else x3line(mm[L3][1-ib],m,mm[L2][1-ia],iall,xyz,dx[L2]);
        if (bug>1) printout("normal","nywt %d nyft %d L %d L2 %d dxL2 %lg %lg %lg\n",
                           nywt,nyft,L,L2,dx[L2][0],dx[L2][1],dx[L2][2]);
        Cross(dx[L],dx[L2],wnorm+3*i);
		}
		else if (nyft==1)  /* cross for normal */
		{
        if (bug>0) printout("normal","warning nywt %d nyft %d not done\n",nywt,nyft);
		}
		else /* nyft==0 an inside edge or corner, look at diagonal points */
		{
        if (bug>0) printout("normal","warning nywt %d nyft %d not done\n",nywt,nyft);
		}
    }
    else if (nywt==2) 
    {
      if (nyft==1) /* cross for normal */
		{ 
        Loop(L,0,3) 
        {
          if (cnear[L][0]=='f') {j=-1; wnear[i]=mm[L][0]; break;}
          else if (cnear[L][1]=='f') {j=1; wnear[i]=mm[L][1]; break;}
        }
        if (j==1) {Cross(dx[(L+1)%3],dx[(L+2)%3],wnorm+3*i);}
        else {Cross(dx[(L+2)%3],dx[(L+1)%3],wnorm+3*i);}
		}
		else if (nyft==2) /* determine what to cross */
		{ 
        Loop(L,0,3) 
        { 
          if (cnear[L][0]=='f' && cnear[L][1]!='w') {j=-1; wnear[i]=mm[L][0]; break;}
          else if (cnear[L][1]=='f' && cnear[L][0]!='w') {j=1; wnear[i]=mm[L][1]; break;}
        }
        if (j==1) {Cross(dx[(L+1)%3],dx[(L+2)%3],wnorm+3*i);}
        else {Cross(dx[(L+2)%3],dx[(L+1)%3],wnorm+3*i);}
		}
		else 
		{
        if (bug>0) printout("normal","warning nywt %d nyft %d not done\n",nywt,nyft);
		}
    }
    else /* insufficient wall segments, look for any near wall point and print warning */
    {
      if (bug>0) printout("normal","warning nywt %d nyft %d not done\n",nywt,nyft);
    }
    wmag=sqrt(wnorm[3*i]*wnorm[3*i]+wnorm[3*i+1]*wnorm[3*i+1]+wnorm[3*i+2]*wnorm[3*i+2]);
    if (wmag>0) Loop(j,0,3) wnorm[3*i+j]/=wmag;
    else /* look for nearby flow points */
    { 
      Loop(L,0,4) inear[L]=ii[L];
      nyft=0;
      if (bug>0) printout("normal"," looking for points near %d %d %d %d\n",ii[0],ii[1],ii[2],ii[3]);
      Loop(ic[0],-1,2) Loop(ic[1],-1,2) Loop(ic[2],-1,2)
      { 
        Loop(L,0,3) inear[L]=max(0,min(i4d[L]-1,ii[L]+ic[L]));
        j=ic[0]+1+3*(ic[1]+1)+9*(ic[2]+1); 
        iflow[j]=In4(inear,i4d); 
        idif[j]=0;
        if (clt[iflow[j]]=='f') { nyft++; Loop(L,0,3) idif[j]+=abs(inear[L]-ii[L]); }
        /*printout("normal","j %d pt %d, %d %d %d clt %c idif %d\n",
         j,iflow[j],inear[0],inear[1],inear[2],clt[iflow[j]],idif[j]); */
      }
      if (nyft==0) notdone++; /* can't find anything */
      else
      { 
        Loop(k,1,4)
        {
          Loop(j,0,27)
          {
            if (idif[j]!=k) continue;
            wnear[i]=iflow[j];
            Loop(L,0,3) wnorm[3*i+L]=xyz[wnear[i]+L*iall]-xyz[m+L*iall];
            wmag=sqrt(wnorm[3*i]*wnorm[3*i]+wnorm[3*i+1]*wnorm[3*i+1]+wnorm[3*i+2]*wnorm[3*i+2]);
            iexpand(wnear[i],i4d,inear);
            if (bug>0) printout("normal"," using %d %d %d %d  wmag %lg \n",inear[0],inear[1],inear[2],inear[3],wmag);
            if (wmag>0) { Loop(L,0,3) wnorm[3*i+L]/=wmag; break;}
          }
          if (wmag>0) break;
        }
      }
    }
    /* if both pt and wnear is on sym plane, make wnorm have no component in sym direction */
    if (nychecksym==1 && wnear[i]>=0)
    { 
      iexpand(whoisw[i],i4d,ii);
      iexpand(wnear[i],i4d,inear);
      Loop(L,0,3) Loop(k,0,2)
      {
        if (ii[L]==k*(i4d[L]-1) && inear[L]==k*(i4d[L]-1) && csym[2*L+k]!='n')
        { 
          m=-1; Loop(n,0,3) if (csym[2*L+k]==cxyz[n]) m=n;
          if (m==-1) continue;
          if (wnorm[3*i+(m+1)%3]==0 && wnorm[3*i+(m+2)%3]==0) continue; /* will get zero if do it */
          if (wnorm[3*i+m]==0) continue; /* don't need to do it */
          wmag=sqrt(wnorm[3*i+(m+1)%3]*wnorm[3*i+(m+1)%3]+wnorm[3*i+(m+2)%3]*wnorm[3*i+(m+2)%3]);
          wnorm[3*i+(m+1)%3]/=wmag;
          wnorm[3*i+(m+2)%3]/=wmag;
          wnorm[3*i+m]=0;
          printout("normal"," sym mod at %d %d %d %d wnorm= %lg %lg %lg\n",ii[0],ii[1],ii[2],ii[3],
                 wnorm[3*i],wnorm[3*i+1],wnorm[3*i+2]);
        }
      }
    }
    if (bug>1) printout("normal","wnear %d wnorm %lg %lg %lg\n",
                       wnear[i],wnorm[3*i],wnorm[3*i+1],wnorm[3*i+2]);
  } /* each independent wall point */
  if (notdone>0) 
  { 
    printout("warning c_wallnorm","WARNING wall norm not set at %d independent wall points\n",notdone);
  }
}

