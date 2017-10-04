/* contains geomcinit, geomcvave, geomcavector, c_geomcprint, geomcgradvol  */
#include "global.h"
/* routine to calculate geometry properties associated with continuity control volumes */

#define Loop(n,a,b) for (n=a;n<b;n++)
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )

static int *i4d,iall;
static double *x[3];
/* --------------------------- */
void geomcinit(void) /* locate arrays for later use */
{ 
  i4d=(int *)need("idim4d");
  iall=Prod4(i4d);
  x[0]=(double *)need("xyz");
  x[1]=x[0]+iall;
  x[2]=x[1]+iall;
}
/* -----------------------------*/
void geomcvave(double *cg, int *ip)    
/*  coefficients for contributions of corners to the volume average */	                         
{
  int i,ig[4],ia[3],in[8],ipt[4];
  int itwo[3]={2,2,2};
  double dx[3][3],vol[8],volb[8],voltot,da[3];
  
  Loop(i,0,8) cg[i]=.125; /* set all equal if can't calculate */
  
  Loop(i,0,3) ig[i]=ip[i]-1;  
  ig[3]=ip[3];   ipt[3]=ip[3];
  /* printout("normal"," ig %d %d %d %d  i4d %d %d %d %d\n",ig[0],ig[1],ig[2],ig[3],i4d[0],i4d[1],i4d[2],i4d[3]); */
  Loop(i,0,3)  if (ig[i]<0 || ig[i] > i4d[i]-2)  return;  /*check legit indices */
  if (ig[3]<0 || ig[3]>i4d[3]-1) return;
  
  Loop3(ia,0,itwo) 
  { 
    Loop(i,0,3) ipt[i]=ig[i]+ia[i];
    in[ia[0]+2*ia[1]+4*ia[2]]=In4(ipt,i4d);
  }
  voltot=0;
  Loop3(ia,0,itwo)
  {
    Loop (i,0,3)
    {
      dx[0][i]=x[i][in[1+2*ia[1]+4*ia[2]]]-x[i][in[2*ia[1]+4*ia[2]]];
      dx[1][i]=x[i][in[ia[0]+2+4*ia[2]]]-x[i][in[ia[0]+4*ia[2]]];
      dx[2][i]=x[i][in[ia[0]+2*ia[1]+4]]-x[i][in[ia[0]+2*ia[1]]];
      /*	 printout("normal"," %lg %lg %lg\n",dx[0][i],dx[1][i],dx[2][i]); */
    }
	 Cross(dx[0],dx[1],da);
	 vol[ia[0]+2*ia[1]+4*ia[2]]=Dot(da,dx[2]);
	 voltot+=  vol[ia[0]+2*ia[1]+4*ia[2]];
  }
  /* Loop(i,0,8) printout("normal"," %lg",vol[i]); printout("normal","\n"); */
  if (voltot<=0) return;
  Loop(i,0,8) vol[i]/=voltot;
  Loop3(ia,0,itwo)
  {
    volb[ia[0]+2*ia[1]+4*ia[2]]=.75*.75*.75*vol[ia[0]+2*ia[1]+4*ia[2]]
    +.75*.75*.25*(vol[1-ia[0]+2*ia[1]+4*ia[2]]+vol[ia[0]+2*(1-ia[1])+4*ia[2]]
                  +vol[ia[0]+2*ia[1]+4*(1-ia[2])])
    +.75*.25*.25*(vol[1-ia[0]+2*(1-ia[1])+4*ia[2]]+vol[ia[0]+2*(1-ia[1])+4*(1-ia[2])]
                  +vol[1-ia[0]+2*ia[1]+4*(1-ia[2])])
    +.25*.25*.25*vol[1-ia[0]+2*(1-ia[1])+4*(1-ia[2])];
  }
  /* Loop(i,0,8) printout("normal"," %lg",volb[i]); printout("normal","\n"); */
  Loop3(ia,0,itwo)
  {
    cg[ia[0]+2*ia[1]+4*ia[2]]=.75*.75*.75*volb[ia[0]+2*ia[1]+4*ia[2]]
    +.75*.75*.25*(volb[1-ia[0]+2*ia[1]+4*ia[2]]+volb[ia[0]+2*(1-ia[1])+4*ia[2]]
                  +volb[ia[0]+2*ia[1]+4*(1-ia[2])])
    +.75*.25*.25*(volb[1-ia[0]+2*(1-ia[1])+4*ia[2]]+volb[ia[0]+2*(1-ia[1])+4*(1-ia[2])]
                  +volb[1-ia[0]+2*ia[1]+4*(1-ia[2])])
    +.25*.25*.25*volb[1-ia[0]+2*(1-ia[1])+4*(1-ia[2])];
  }
  /* Loop(i,0,8) printout("normal"," %lg",cg[i]); printout("normal","\n"); */
}
/* ------------------------------ */
/*  for a cont. surface, portion of area vector associated with each corner point */  
/* ic is lower corner grid index */
void geomcavector(double (*area)[3], int *ic, int L)  
{
  int ig[4],ipt[2][2],ia,ib,L2,L3,i,j;
  double dxa[2][3],dxb[2][3],a[2][2][3],ab[2][2][3];
  
  L2=(L+1)%3, L3=(L+2)%3;
  ig[3]=ic[3]; ig[L]=ic[L];
  Loop(i,0,3) Loop(j,0,4) area[j][i]=0;
  
  Loop(ia,0,2) Loop(ib,0,2)
  { 
    ig[L2]=ic[L2]+ia; ig[L3]=ic[L3]+ib;
	 ipt[ia][ib]=In4(ig,i4d);
	 if (ipt[ia][ib]<0 || ipt[ia][ib]>iall) return;   /* out of range indices */
	 /*	 printout("normal"," ia %d ib %d ig %d %d %d ipt %d %lg %lg %lg\n",
     ia,ib,ig[0],ig[1],ig[2],ipt[ia][ib],x[0][ipt[ia][ib]],x[1][ipt[ia][ib]],x[2][ipt[ia][ib]]); */
  }
  Loop(j,0,2)
  {
    Loop(i,0,3)
    {
      dxa[j][i]=.5*(x[i][ipt[1][j]]-x[i][ipt[0][j]]);
		dxb[j][i]=.5*(x[i][ipt[j][1]]-x[i][ipt[j][0]]);
    }
    /*	printout("normal"," j %d dxa %lg %lg %lg   dxb %lg %lg %lg \n",
     j,dxa[j][0],dxa[j][1],dxa[j][2],dxb[j][0],dxb[j][1],dxb[j][2]); */
  }
  Loop(ia,0,2) Loop(ib,0,2) 
  {
    Cross(dxa[ib],dxb[ia],a[ia][ib]);    /* quarter area based on corners */
	 /*	 printout("normal"," ia ib %d %d area %lg %lg %lg\n",ia,ib,a[ia][ib][0],a[ia][ib][1],a[ia][ib][2]); */
  }
  Loop(ia,0,2) Loop(ib,0,2) 
  {
    Loop(i,0,3)  /* real 'quarter' area */
    ab[ia][ib][i]=.75*.75*a[ia][ib][i]+.75*.25*(a[1-ia][ib][i]+a[ia][1-ib][i])
    +.25*.25*a[1-ia][1-ib][i];
	 /*	 printout("normal"," ia ib %d %d areab %lg %lg %lg\n",ia,ib,ab[ia][ib][0],ab[ia][ib][1],ab[ia][ib][2]); */
  }
  Loop(ia,0,2) Loop(ib,0,2) Loop(i,0,3)   /* corner point interpolation area */
  area[ia+2*ib][i]=.75*.75*ab[ia][ib][i]+.75*.25*(ab[1-ia][ib][i]+ab[ia][1-ib][i])
  +.25*.25*ab[1-ia][1-ib][i];
}
/* --------------------- */
/* command level print checks */
void c_geomcprint(FILE *fpin, FILE *fprint)
{ 
  int i,ip[4],is[4],ie[4],ig[4],L;
  double cg[8],cgsum,area[4][3],vol,grad[8][3];
  
  geomcinit();
  Loop(i,0,4) 
  { 
    is[i]=readint(fpin);
	 ie[i]=readint(fpin)+1;
  }
  Loop(ig[3],is[3],ie[3]) Loop(ig[2],is[2],ie[2])
  Loop(ig[1],is[1],ie[1]) Loop(ig[0],is[0],ie[0])
  { 
    Loop(i,0,3) ip[i]=ig[i]+1; ip[3]=ig[3];
	 geomcvave(cg,ip);
	 printout("print","%d %d %d %d \n   cg ",ig[0],ig[1],ig[2],ig[3]);
	 cgsum=0;
	 Loop(i,0,8) {printout("print"," %lg",cg[i]); cgsum+=cg[i]; }
	 printout("print"," cg-sum  %lg\n",cgsum);
	 Loop(L,0,3)
    {
      geomcavector(area,ig,L);
      printout("print"," L %d areas",L);
      Loop(i,0,4) printout("print","    %lg %lg %lg",area[i][0],area[i][1],area[i][2]);
      printout("print","\n");
    }
	 vol=geomcgradvol(ip,grad);
	 printout("print","gradvol, vol=%lg\n",vol);
	 Loop(L,0,3) 
    {
      printout("print","d/dx%d",L);
      Loop(i,0,8) printout("print"," %lg",grad[i][L]);
      printout("print","\n");
    }
  }			 
}
/* ------------------- */
/* returns cont. c.v. volume, and grad[8][3] eight corners,i,j average for each surface */
/* grad*vol = di x dj d/dk + dj x dk d/di + dk x di d/dj. vol=dk dot (di x dj) */
/* sum for 8 sections to get some non-parallelogram influence */
/* ip is the pressure index associated with the cont c.v. */
double geomcgradvol(int *ip, double (*grad)[3])
{
  int i,j,k,L,L2,L3,ia[3],ib[3],ig[4];
  int two[3]={2,2,2};
  double xc[8][3],fmid[3][2],dd[8][3],da[3][3];
  double vol,dx[3][3]; /* [i,j,or k directions][x,y,z components] */
  
  /* printout("normal","ip %d %d %d %d\n",ip[0],ip[1],ip[2],ip[3]); */
  Loop(i,0,8) Loop(j,0,3) grad[i][j]=0; 
  vol=0;
  Loop(i,0,3) if (ip[i]==0 || ip[i]>i4d[i]-1) return vol; /* out or range */
  
  ig[3]=ip[3]; /* x at 8 corners  */
  Loop3(ia,0,two) 
  {
    Loop(i,0,3) ig[i]=ip[i]-1+ia[i];
    j=ia[0]+2*ia[1]+4*ia[2];
    k=In4(ig,i4d); 
    Loop(i,0,3) xc[j][i]=x[i][k]; 
    /* printout("normal","ia %d %d %d xc %lg %lg %lg\n",ia[0],ia[1],ia[2],xc[j][0],xc[j][1],xc[j][2]); */
  }
  Loop3(ia,0,two)
  {
    Loop(i,0,3) {fmid[i][1]=.25+.5*ia[i]; fmid[i][0]=1-fmid[i][1]; }
    Loop(i,0,8) Loop(j,0,3) dd[i][j]=0;
    Loop(L,0,3) 
    { 
      L2=(L+1)%3; L3=(L+2)%3;
      Loop(ib[L2],0,2) Loop(ib[L3],0,2)
	   { 
        ib[L]=1;
	     j=ib[0]+2*ib[1]+4*ib[2];
	     dd[j][L]+=fmid[L2][ib[L2]]*fmid[L3][ib[L3]];
	     ib[L]=0; 
	     j=ib[0]+2*ib[1]+4*ib[2];
	     dd[j][L]-=fmid[L2][ib[L2]]*fmid[L3][ib[L3]];
	   }
    }
    Loop(i,0,3) Loop(j,0,3) 
    {
      dx[i][j]=0;
      Loop(k,0,8) dx[i][j]+=dd[k][i]*xc[k][j];
    }
    Cross(dx[0],dx[1],da[2]);
    Cross(dx[1],dx[2],da[0]);
    Cross(dx[2],dx[0],da[1]);
    vol+=Dot(dx[0],da[0]);
    Loop3(ib,0,two)
    {
      k=ib[0]+2*ib[1]+4*ib[2];
      Loop(i,0,3) Loop(j,0,3) grad[k][j]+= dd[k][i]*da[i][j];
    }
  }
  if (vol>0) Loop(k,0,8) Loop(j,0,3) grad[k][j]/=vol;
  return vol/8.;
}
