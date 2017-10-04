/* contains whoelseadd, c_gridmatch */
#include "global.h"

#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Equal(a,b) ( abs((a)-(b))<=tol*max(abs(a),abs(b)) )
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Matchi(j) (mat[0][j][0]==mat[1][j][0] &&  mat[0][j][1]==mat[1][j][1] )
#define Matchj(i) ( mat[i][0][0]==mat [i][1][0] && mat[i][0][1]==mat[i][1][1] )

/* -------------------------*/
void whoelseadd(int **whoelse, int k, int i)
{	
  int *whoelsett,j;
  if (whoelse[k]==0)
  {
    whoelse[k]=(int *)smalloca(2,'i');
    whoelse[k][0]=1;
    whoelse[k][1]=i;
  }
  else
  { 
    whoelsett=(int *)smalloca(whoelse[k][0]+2,'i');
    Loop(j,0,whoelse[k][0]) whoelsett[j+1]=whoelse[k][j+1];
    whoelsett[0]=whoelse[k][0]+1;
    whoelsett[whoelse[k][0]+1]=i;
    free(whoelse[k]);
    whoelse[k]=whoelsett;
  }
}
/* ----------------------- */
/* set integer arrays  match, whoelse, idim4dp, matchpc
     and character array matchside(iallp,6) 
*/
void c_gridmatch(FILE *fpin, FILE *fprint)
{  
  int *i4d; double *xyz;  char *csym; /* needed arrays */
  double *roundoff; /* use if available */
  int *match, *nogpts,**whoelse;   /* arrays set */
  int *matchpc,*noppts, *i4dp; 
  char *cmatch;
  
  double *x,*y,*z, tm=0,tn=0;
  int iall,i,j=0,k,n,m,ii[4],jj[4],L;
  
  int icorrect,matchmin,nn[3];
  int iallp,ip[4],ia[4],ma,mb,mc,na[4],nb[4],mat[2][2][2],nyok[4];
  int two[3]={2,2,2};
  int L2,L3,L4,ncheck,ncheckold;
  double tol=1.e-11;
  
  int nyprint=0; /*  to debug match */
  int bug=0; /* to debug matchpc */
  
  /* nyprint=readint(fpin);  */
  
  /* set match */
  i4d=(int *)need("idim4d");
  csym=(char *)need("csym");
  iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
  roundoff=(double *)find("roundoff");
  if (roundoff>0) tol=roundoff[0];	
  xyz=(double *)need("xyz");
  x=xyz;
  y=xyz+iall;
  z=xyz+2*iall;
  
  nogpts=(int*)createarray("nogpts",1,'i',0);
  nogpts[0]=iall;
  match=(int *)createarray("match",iall,'i',0);
  Loop(i,0,iall) match[i]=i;
  whoelse=(int **)createarray("whoelse",iall,'p',0);
  Loop(i,0,iall) whoelse[i]=0;
  
  i4dp=(int *)createarray("idim4dp",4,'i',0);
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  iallp=Prod4(i4dp);
  noppts=(int *)createarray("noppts",1,'i',0);
  noppts[0]=iallp;
  cmatch=(char *)createarray("matchside",iallp*6,'c',0);
  Loop(i,0,iallp*6) cmatch[i]='n';
  matchpc=(int *)createarray("matchpc",iallp*2,'i',0);
  Loop(i,0,2*iallp) matchpc[i]=-1;   /* default not sleeping */
  
  icorrect=0;
  /* check along i-lines,  j-lines, and k lines  */
  Loop(ii[3],0,i4d[3]) Loop(ii[2],0,i4d[2]) Loop(ii[1],0,i4d[1]) Loop(ii[0],0,i4d[0])
  {
    m=In4(ii,i4d);
    tm=ii[3];
    Loop(i,0,4) jj[i]=ii[i];
    matchmin=m; Loop(L,0,3) nn[L]=-1;
    Loop(L,0,3) 
    { 
      if (jj[L]==0) continue; 
      jj[L]--; n=In4(jj,i4d);
      tn=jj[3];
      jj[L]++;
      if (nyprint>2) printout("normal","ii %d %d %d %d m %d n %d xyz %lg %lg    %lg %lg    %lg %lg    %lg %lg\n",
                             ii[0],ii[1],ii[2],ii[3],m,n,x[m],x[n],y[m],y[n],z[m],z[n],tm,tn);
      
      if (x[m]==x[n] && y[m]==y[n] && z[m]==z[n]&& tm==tn) nn[L]=n;
      
      else if (Equal(x[m],x[n]) && Equal(y[m],y[n]) && Equal(z[m],z[n]) && Equal(tm,tn) )
      { 
        nn[L]=n;
        if (nyprint>2) printout("normal","fixroundoff m xyz %d %lg %lg %lg n xyz %d %lg %lg %lg\n",
                               m,x[m],y[m],z[m],n,x[n],y[n],z[n]);
        x[m]=x[n]; y[m]=y[n]; z[m]=z[n]; icorrect++;
      }
    }
    Loop(L,0,3) 
    if (nn[L]>=0) 
    { 
      while (match[nn[L]]!=nn[L]) nn[L]=match[nn[L]];
      matchmin=min(matchmin,nn[L]);
    }
    if (matchmin<m) 
    { 
      match[m]=matchmin;
      Loop(L,0,3) if (nn[L]>=0) match[nn[L]]=matchmin;
      if (nyprint>2) printout("normal"," matchmin %d for %d %d %d %d\n",matchmin,m,nn[0],nn[1],nn[2]);
    }
  }
  if (icorrect>0) printout("normal"," roundoff matches corrected at %d points\n",icorrect);
  
  
  /* match repeat information, set 2-d repeat only here  full repeat set later */
  Loop(i,0,i4d[0]) Loop(j,0,i4d[1]) Loop(k,0,i4d[2]) Loop(n,0,i4d[3])
  { 
    m=i+i4d[0]*(j+i4d[1]*(k+i4d[2]*n));
    if (csym[6]=='R' && i==i4d[0]-1) match[m]=match[i4d[0]*(j+i4d[1]*(k+i4d[2]*n))];
    if (csym[8]=='R' && j==i4d[1]-1) match[m]=match[i+i4d[0]*(i4d[1]*(k+i4d[2]*n))];
    if (csym[10]=='R' && k==i4d[2]-1) match[m]=match[i+i4d[0]*(j+i4d[1]*(i4d[2]*n))];
    /* if (nyrepeat[3]==1 && n==i4d[3]-1) match[m]=match[i+i4d[0]*(j+i4d[1]*(k))]; time repeat omitted */
  }
  /* consolidate match */
  Loop (j,0,10)
  { 
    icorrect=0;
    Loop(n,0,iall)
    { 
      m=match[n];
      if (m==n) continue;
      if (match[m]!=m) 
      { 
        if (nyprint>1) printout("normal"," n match[n] match[match[n]] %d %d %d\n",
                               n,m,match[m]);
        match[n]=min(m,match[m]); 
        match[m]=match[n]; icorrect++; 			
      }		
    }
    printout("normal"," match consolidation for %d points\n",icorrect);
    if (icorrect==0) break;
  }
  /* initial setting of whoelse   */
  Loop(i,0,iall) 
  if (match[i]!=i) whoelseadd(whoelse,match[i],i);
  
  /* set matchpc (repeat bndry on match set after matchpc) */
  /**********************/
  
  /* check for sleeping internal volumes and set  matchpc  i,j,k and ij corners for now*/
  ncheck=0;
  Loop(ip[3],0,i4d[3]) 
  { 
    ii[3]=ip[3],na[3]=ip[3],nb[3]=ip[3];
    /* start with i lines j lines and k lines */
    Loop(ip[2],0,i4dp[2]) Loop(ip[1],0,i4dp[1]) Loop(ip[0],0,i4dp[0])
    { 
      if (bug>1) printout("normal","ip %d %d %d %d \n",ip[0],ip[1],ip[2],ip[3]);
      if (bug>1) printout("normal","  mpc %d\n",matchpc[In4(ip,i4dp)]);
      if (matchpc[In4(ip,i4dp)]<0)
      { 
        Loop(i,0,3) {na[i]=ip[i]; nb[i]=ip[i]; }
        if (bug>1) printout("normal","ip %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]);
        Loop3(ia,0,two)    
        {
          Loop(i,0,3) ii[i]=max(0,min(i4d[i]-1,ip[i]-1+ia[i]));
          mat[ia[0]][ia[1]][ia[2]]=match[In4(ii,i4d)];
        }
        if (bug>1) printout("normal"," mat set\n");
        if (csym[6]!='R' && ip[0]>0 && ip[0]<i4d[0])  /* check i matches */
          while (Matchi(0) && Matchi(1) && nb[0]<i4d[0])
          { 
            nb[0]++; ia[0]=1; ncheck++;
            Loop(ia[1],0,2) Loop(ia[2],0,2)
            {
              Loop(i,0,3) ii[i]=max(0,min(i4d[i]-1,nb[i]-1+ia[i]));
              mat[1][ia[1]][ia[2]]=match[In4(ii,i4d)];
            }
          }
        if (nb[0]>na[0])  /* have i matches */
        {
          na[0]--; ma=In4(na,i4dp); mb=In4(nb,i4dp);
          i=na[0]+1; j=nb[0];
          if (bug>0) printout("normal","imatch ma mb %d %d i %d %d j %d k %d\n",ma,mb,i,j,na[1],na[2]);
          Loop(nb[0],i,j)
          {
            mc=In4(nb,i4dp);
            matchpc[mc]=ma;
            matchpc[mc+iallp]=mb;
            cmatch[mc]='I'; cmatch[mc+iallp]='i';
          }
          continue;
        }
        if (bug>1) printout("normal"," i %d %d %d %d\n",i4d[0],i4d[1],i4d[2],i4d[3]);
        if (csym[8]!='R' && ip[1]>0 && ip[1]<i4d[1])  /* check j matches */
          while (Matchj(0) && Matchj(1) && nb[1]<i4d[1])
          { 
            nb[1]++; ia[1]=1; ncheck++;
            Loop(ia[0],0,2) Loop(ia[2],0,2)
            {
              Loop(i,0,3) ii[i]=max(0,min(i4d[i]-1,nb[i]-1+ia[i]));
              mat[ia[0]][1][ia[2]]=match[In4(ii,i4d)];
            }
          }
        if (nb[1]>na[1])  /* have j matches */
        { 
          na[1]--; ma=In4(na,i4dp); mb=In4(nb,i4dp);
          i=na[1]+1; j=nb[1];
          if (bug>0) printout("normal","jmatch ma mb %d %d i %d j %d %d k %d\n",ma,mb,na[0],i,j,na[2]);
          Loop(nb[1],i,j)
          {
            mc=In4(nb,i4dp);
            matchpc[mc]=ma;
            matchpc[mc+iallp]=mb;
            cmatch[mc+2*iallp]='J'; cmatch[mc+3*iallp]='j';
          }
          continue;
        }
        if (bug>1) printout("normal"," j\n");
        if (csym[10]!='R' && ip[2]>0 && ip[2]<i4d[2])  /* check k matches */
          while (nb[2]<i4d[2] && mat[0][0][0]==mat[0][0][1] && mat[1][0][0]==mat[1][0][1]
                 && mat[0][1][0]==mat[0][1][1] && mat[1][1][0]==mat[1][1][1])
          { nb[2]++; ia[2]=1; ncheck++;
            Loop(ia[0],0,2) Loop(ia[1],0,2)
            {Loop(i,0,3) ii[i]=max(0,min(i4d[i]-1,nb[i]-1+ia[i]));
              mat[ia[0]][ia[1]][1]=match[In4(ii,i4d)];
            }
          }
        if (nb[2]>na[2])  /* have k matches */
        { 
          na[2]--; ma=In4(na,i4dp); mb=In4(nb,i4dp);
          i=na[2]+1; j=nb[2];
          if (bug>0) printout("normal","kmatch ma mb %d %d i %d j %d k %d %d\n",ma,mb,na[0],na[1],i,j);
          Loop(nb[2],i,j)
          {
            mc=In4(nb,i4dp);
            matchpc[mc]=ma;
            matchpc[mc+iallp]=mb;
            cmatch[mc+4*iallp]='K'; cmatch[mc+5*iallp]='k';
          }
          continue;
        }
        if (bug>1) printout("normal"," k\n ");
      }
    }
    if (bug>0) printout("normal","done i j k matches\n");
    /* now look for i-j corners */
    Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
    {
      mc=In4(ip,i4dp);
      if (matchpc[mc]<0)
      { 
        Loop(i,0,3) {na[i]=ip[i]; nb[i]=ip[i]; }
        Loop3(ia,0,two)    
        {
          Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
          mat[ia[0]][ia[1]][ia[2]]=match[In4(ii,i4d)];
        }
        if (Matchi(0) && Matchj(0)) 
        {na[0]++; nb[1]++; cmatch[mc+iallp]='J'; cmatch[mc+3*iallp]='I';  }
        else if (Matchi(1) && Matchj(0)) 
        {na[1]--; nb[0]++; cmatch[mc+iallp]='j'; cmatch[mc+2*iallp]='I';}
        else if (Matchi(0) && Matchj(1)) 
        { na[0]--; nb[1]++; cmatch[mc]='J'; cmatch[mc+3*iallp]='i';}
        else if (Matchi(1) && Matchj(1)) 
        {na[1]--; nb[0]--; cmatch[mc]='j'; cmatch[mc+2*iallp]='i';}
        else continue;
        ma=In4(na,i4dp);
        mb=In4(nb,i4dp);
        /* printout("normal","corner ip %d %d na %d %d nb %d %d matchpcip %d %d  cmatch",
         ip[0],ip[1],na[0],na[1],nb[0],nb[1],matchpc[mc],matchpc[mc+iallp]);
         Loop(i,0,6) printout("normal"," %c",cmatch[mc+i*iallp]); printout("normal","\n");
         printout("normal","ma %d matchpc %d %d ",ma, matchpc[ma],matchpc[ma+iallp]); 
         printout("normal","  mb %d matchpc %d %d ",mb, matchpc[mb],matchpc[mb+iallp]);  */
        matchpc[mc]=ma; 
        matchpc[mc+iallp]=mb;
        if (matchpc[ma]<0 || matchpc[mb]<0) continue;
        if (matchpc[ma]==matchpc[mb])
        { 
          matchpc[mc]=matchpc[ma+iallp]; 
          matchpc[mc+iallp]=matchpc[mb+iallp];
        }
        else if (matchpc[ma+iallp]==matchpc[mb])
        {
          matchpc[mc]=matchpc[ma]; 
          matchpc[mc+iallp]=matchpc[mb+iallp];
        }
        else if (matchpc[ma]==matchpc[mb+iallp])
        { 
          matchpc[mc]=matchpc[ma+iallp]; 
          matchpc[mc+iallp]=matchpc[mb];
        }
        else if (matchpc[ma+iallp]==matchpc[mb+iallp])
        { 
          matchpc[mc]=matchpc[ma]; 
          matchpc[mc+iallp]=matchpc[mb];
        }
        ncheck++;
        /* printout("normal","   mc %d matchpc %d %d\n",mc,matchpc[mc],matchpc[mc+iallp]); */
      }
    }
  }  /* ip[3] */
  
  printout("normal","%d sleeping internal continuity volumes\n",ncheck);
  
  /* pressure boundaries */
  Loop(L,0,3) 
  {
    L2=(L+1)%3; L3=(L+2)%3;
    Loop(na[L2],0,i4d[L2]+1) Loop(na[L3],0,i4d[L3]+1) Loop(na[3],0,i4d[3])
    {  
      if (csym[6+2*L]=='r' || csym[6+2*L]=='R')
      {
        na[L]=1; ma=In4(na,i4dp);
        na[L]=i4d[L]-1; mb=In4(na,i4dp);
        na[L]=0; mc=In4(na,i4dp);
        matchpc[mc]=ma; matchpc[mc+iallp]=mb;
        cmatch[mc+(1+2*L)*iallp]='r';
        na[L]=i4d[L]; mc=In4(na,i4dp);
        matchpc[mc]=mb; matchpc[mc+iallp]=ma;
        cmatch[mc+2*L*iallp]='R';
        continue;
      }
      if (csym[6+2*L]=='s')
      { 
        na[L]=1; ma=In4(na,i4dp);
        na[L]=0; mc=In4(na,i4dp);
        matchpc[mc]=ma; matchpc[mc+iallp]=ma;
      }
      if (csym[7+2*L]=='s')
      { 
        na[L]=i4d[L]-1; mb=In4(na,i4dp);
        na[L]=i4d[L]; mc=In4(na,i4dp);
        matchpc[mc]=mb; matchpc[mc+iallp]=mb;
      }
    }
  }
  
  /* consolidate matchpc, allow only in currently differing direction */
  ncheckold=0;
  Loop(k,0,10)
  { 
    ncheck=0; 
    Loop(ip[3],0,i4d[3]) Loop(ip[2],0,i4dp[2]) Loop(ip[1],0,i4dp[1]) Loop(ip[0],0,i4dp[0])
    {
      i=In4(ip,i4dp);
      if (bug>1) printout("normal","matchpc at %d %d %d %d,  %d\n", ip[0],ip[1],ip[2],ip[3],matchpc[i]);
      if (matchpc[i]<0) continue; 
      /* omit self */
      if (matchpc[i]==i) matchpc[i]=matchpc[i+iallp];
      if (matchpc[i+iallp]==i) matchpc[i+iallp]=matchpc[i];
      iexpand(matchpc[i],i4dp,na); 
      iexpand(matchpc[i+iallp],i4dp,nb);
      Loop(L,0,4) {nyok[L]=0; if (ip[L]!=na[L] || ip[L]!=nb[L]) nyok[L]=1; }
      if (bug>1) 
      { 
        iexpand(matchpc[i],i4dp,na); iexpand(matchpc[i+iallp],i4dp,nb);
        printout("normal","                 from %d %d %d %d and %d %d %d %d\n",
                na[0],na[1],na[2],na[3],nb[0],nb[1],nb[2],nb[3]);
      }
      /* consolicate */
      j=matchpc[i];
      if (matchpc[j]>=0) 
      { 
        iexpand(matchpc[j],i4dp,na);
        n=1; 
        Loop(L,0,4) {if (ip[L]!=na[L] && nyok[L]==0) {n=0; break;}}
        if (n==1)
        {
          matchpc[i]=matchpc[j]; 
          if (matchpc[i+iallp]==matchpc[i]) matchpc[i+iallp]=matchpc[j+iallp];
          if (bug>1) 
          { 
            iexpand(matchpc[i+iallp],i4dp,nb);
            printout("normal","            to  %d %d %d %d and %d %d %d %d\n",
                    na[0],na[1],na[2],na[3],nb[0],nb[1],nb[2],nb[3]);
          }
          ncheck++;
        }
      }
      j=matchpc[i+iallp];
      if (matchpc[j]>=0) 
      {  
        iexpand(matchpc[j],i4dp,na);
        n=1; 
        Loop(L,0,4) {if (ip[L]!=na[L] && nyok[L]==0) {n=0; break;}}
        if (n==1)
        {
          matchpc[i+iallp]=matchpc[j+iallp]; 
          if (matchpc[i+iallp]==matchpc[i]) matchpc[i]=matchpc[j];
          if (bug>1) 
          { 
            iexpand(matchpc[i],i4dp,na); iexpand(matchpc[i+iallp],i4dp,nb);
            printout("normal","            to  %d %d %d %d and %d %d %d %d\n",
                    na[0],na[1],na[2],na[3],nb[0],nb[1],nb[2],nb[3]);
          }
          ncheck++;}
      }
      /* check */
		/*	if (k==9 || ncheck==0)
       if (matchpc[matchpc[i]]>=0 || matchpc[matchpc[i+iallp]]>=0)
       { printout("normal","warning: at %d %d %d %d, ip %d matchpc %d %d",
       ip[0],ip[1],ip[2],ip[3],i,matchpc[i],matchpc[i+iallp]);
       printout("normal","matchpc of those %d %d,  %d %d\n",
       matchpc[matchpc[i]],matchpc[matchpc[i]+iallp],
       matchpc[matchpc[i+iallp]],matchpc[matchpc[i+iallp]+iallp]);
       }
		 */
    } /* ip loop */
    printout("normal","matchpc simplifcation at %d points\n",ncheck);
    if (ncheck==0) break;
    if (ncheck==ncheckold) k=max(k,8);		
    ncheckold=ncheck;
  }
  if (ncheck>0) 
  { 
    printout("error c_gridmatch","matchpc simplification failed, probable grid error, exit\n");
    exitm4d(0);
  }
  
  /************************/
  /* match repeat information */
  
  Loop(L,0,3)   /* Loop(L,0,4) omit time repeat */
  if (csym[6+2*L]=='r')
  { 
    L2=(L+1)%4; L3=(L+2)%4; L4=(L+3)%4;
    Loop(ii[L2],0,i4d[L2]) Loop(ii[L3],0,i4d[L3]) Loop(ii[L4],0,i4d[L4]) 
    { 
      ii[L]=i4d[L]-1; m=In4(ii,i4d);
      ii[L]=0; n=In4(ii,i4d);
      if (match[m]!=match[n])
      {	
        whoelseadd(whoelse,match[n],match[m]);
        if (whoelse[match[m]]>0) 
        { 
          Loop(j,0,whoelse[match[m]][0]) 
          {
            whoelseadd(whoelse,match[n],whoelse[match[m]][j+1]);
          }
          free (whoelse[match[m]]);
          whoelse[match[m]]=0;
        }
      }
      Loop(j,0,whoelse[match[n]][0]) match[whoelse[match[n]][j+1]]=match[n];
    }
    printout("normal","match repeating boundary set for L= %d\n",L);
  }
  if (nyprint>0) /* print whoelse */
  { 
    k=0;
    Loop(i,0,iall)
    if (whoelse[i]>0)
    { 
      k++;
      printout("normal","\n mult points at %d  ",i);
      Loop(j,0,whoelse[i][0]) 
      printout("normal"," %d",whoelse[i][j+1]);		
    }
    printout("normal","\n");
    
    printout("normal"," %d multiple points\n",k);
  }
}
