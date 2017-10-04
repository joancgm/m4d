/* contains c_gradprop */
#include "global.h"
#include "psleepcalc.h"

#define Loop(n,a,b) for (n=a;n<b;n++) 
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])

/* set a mid-point gradient arrays from on the points arrays, up to 20 names */
/* use roundoff filter if roundoffd exists */
/*------------------ -------------------*/
void c_gradprop(FILE *fpin, FILE *fprint)
{  
  int *i4d,*wherep; /* need */
  double *cpsleepm=0,*roundoff=0;  /* use if exists */
  char *name,*nameto[3],ndi[5]="ddxi",ndd[3][5]={"ddx1","ddx2","ddx3"};
  double *vf,*vto[3],*vpto=0,vround;
  
  int i,j,k,n,nper,size,i4dm[4],i4dp[4],iall,iallm,iallp,ip[4],mmid,ii[4],mc[8],ia[4];
  int two[3]={2,2,2};
  double volc,grad[8][3];
  
  i4d=(int *)need("idim4d");
  wherep=(int *)need("wherep");
  Loop(i,0,3) i4dp[i]=i4d[i]+1; i4dp[3]=i4d[3];
  Loop(i,0,3) i4dm[i]=i4d[i]-1; i4dm[3]=i4d[3];
  iall=Prod4(i4d);
  iallp=Prod4(i4dp);
  iallm=Prod4(i4dm);
  cpsleepm=(double *)find("cpsleepm");
  roundoff=(double *)find("roundoffd");
  if (cpsleepm>0) vpto=(double *)tmalloca(iallp,'d');
  geomcinit();
  
  nper=readint(fpin);  /* =1 for all in grad, ddxi,  =3 for ddx1,ddx2,ddx3 arrays */
  if (nper==1) printout("normal"," creating arrays ___%s\n",ndi);
  else if (nper==3) printout("normal"," creating arrays ___%s ___%s ___%s\n",ndd[0],ndd[1],ndd[2]);
  else 
  {
    printout("error c_gradprop"," input error, nper=%d, needs to be 1 or 3\n",nper); 
    exitm4d(0); 
  }
  
  while(1)
  { 
    name=readname(fpin);
    if (name[0]=='\0')   break;
    if (strcmp(name,"c:")==0) { fseek(fpin,(long)(-3),1);  break; }
    size=arraysize(name);
    if (size!=iall) 
    { 
      printout("warning c_gradprop"," input error array %s, size %d not an on the points array\n ",name,size); break;}
    vf=(double *)need(name);
    Loop(i,0,nper) 
    { 
      nameto[i]=(char *)tmalloca(strlen(name)+5,'c'); 
      Loop(j,0,(int)strlen(name)) nameto[i][j]=name[j];
    }
    if (nper==1) 
    { 
      Loop(n,0,5) nameto[0][strlen(name)+n]=ndi[n]; 
      nameto[0][strlen(name)+5]='\0';
      printout("normal","from %s creating %s\n",name,nameto[0]);
      vto[0]=(double *)createarray(nameto[0],iallm*3,'d',0);
      Loop(n,0,iallm*3) vto[0][n]=0;
      vto[1]=vto[0]+iallm; vto[2]=vto[1]+iallm;
    }
    else 
    {  
      Loop(j,0,3) 
      {
        Loop(n,0,5) nameto[j][strlen(name)+n]=ndd[j][n]; 
        nameto[j][strlen(name)+5]='\0';
        vto[j]=(double *)createarray(nameto[j],iallm,'d',0);
        Loop(n,0,iallm) vto[j][n]=0;
      }
      printout("normal","from %s creating %s %s and %s\n",name,nameto[0],nameto[1],nameto[2]);
    }
    /* loop over contiuity control volumes */
    Loop(ip[3],0,i4d[3]) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
    if (wherep[In4(ip,i4dp)]>=0)   /* do only for non-zero cont c.v. */
    { 
      mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3])); 	
      volc=geomcgradvol(ip,grad);   /* geometry arrays */
      ii[3]=ip[3];        /* corner indices */
      Loop3(ia,0,two)
      { 
        Loop(i,0,3) ii[i]=ip[i]-1+ia[i];
        mc[ia[0]+2*ia[1]+4*ia[2]]=In4(ii,i4d);
      }         
      Loop(j,0,3)   /*   gradient */
      {
        Loop(k,0,8) vto[j][mmid] += grad[k][j]*vf[mc[k]];
        if (roundoff>0)
        {
          vround=0;
          Loop(k,0,8) vround+=abs(grad[k][j])*vf[mc[k]];
          if (abs(vto[j][mmid])<roundoff[0]*abs(vround)) vto[j][mmid]=0;                        
        }
      }
    }
    if (cpsleepm>0) /* put on p grid, calc sleep and copy back */
      Loop(j,0,3)
    {
      Loop(i,0,iallp) vpto[i]=0;
      Loop(ip[3],0,i4d[3]) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
      { 
        mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3])); 	
        i=In4(ip,i4dp);
        vpto[i]=vto[j][mmid];
      }
      psleepcalc(fprint,vpto,cpsleepm);
      Loop(ip[3],0,i4d[3]) Loop(ip[2],1,i4d[2]) Loop(ip[1],1,i4d[1]) Loop(ip[0],1,i4d[0])
      { 
        mmid=ip[0]-1+i4dm[0]*(ip[1]-1+i4dm[1]*(ip[2]-1+i4dm[2]*ip[3])); 	
        i=In4(ip,i4dp);
        vto[j][mmid]=vpto[i];
      }
    }
  }
}