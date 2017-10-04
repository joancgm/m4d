/* contains c_mirror */
/* start with coding mirror option, note   character type not tested!!*/
#include "global.h"

#define Loop(n,a,b) for (n=a;n<b;n++)
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))

void c_mirror(fpin,fprint)
FILE *fpin,*fprint;
{ char todo, aborc,*append, *namefrom,*nameto;
	int side;
	
	int *i4d;
	int i,j,k,idnew[4],idold[4],L,iall,iallp,iallm,iadd,irep,jrep,sign,ii[4],iiold[4],L2,L3,iallold,iallnew;
	Array *af;
	double *v,*vf,vv;  char *c,*cf,cc;
	
	i4d=(int *)need("idim4d");
	iall=i4d[0]*i4d[1]*i4d[2]*i4d[3];
	iallp=(i4d[0]+1)*(i4d[1]+1)*(i4d[2]+1)*i4d[3];
	iallm=(i4d[0]-1)*(i4d[1]-1)*(i4d[2]-1)*i4d[3];
	append=readname(fpin);
	todo=read1charname(fpin);
	if (todo=='m') /* mirror */
	{  aborc=read1charname(fpin);
		L=-1; if (aborc=='a') L=0; else if (aborc=='b') L=1; else if (aborc=='c') L=2;
		L2=(L+1)%3; L3=(L+2)%3;
		side=readint(fpin);
		printout("normal","append %s to mirror %c (L=%d) to side %d\n",append,aborc,L,side);
		if (L<0) 
      {
        printout("error c_mirror"," input error, exit\n"); 
        exitm4d(0);
      }
		while (1)
		{ namefrom=readname(fpin);
			af=findarray(namefrom);
			if (af==0) {printout("normal"," %s ending regrid\n",namefrom); break; }
			else printout("normal"," from %s",namefrom);
			if (af->size==iall) {iadd=0; irep=1;}
			else if (af->size==iallp) {iadd=1; irep=1;}
			else if (af->size==iallm) {iadd=-1; irep=1;}
			else if (af->size%iall==0) {iadd=0; irep=af->size/iall;}
			else if (af->size%iallp==0) {iadd=1; irep=af->size/iallp;}
			else if (af->size%iallm==0) {iadd=-1; irep=af->size/iallm;}
			else 
			{printout("warning c_mirror"," %s size %d is not standard, cannot mirror\n",namefrom,af->size);
				continue;
			}
			Loop(i,0,3) idold[i]=i4d[i]+iadd; idold[3]=i4d[3];
			iallold=idold[0]*idold[1]*idold[2]*idold[3];
			Loop(i,0,4) idnew[i]=idold[i];
			if (iadd==0) idnew[L]+=i4d[L]-1;
			else if (iadd==1) idnew[L]=2*(idnew[L]-1);
			else if (iadd==-1) idnew[L]=2*idnew[L];
			iallnew=idnew[0]*idnew[1]*idnew[2]*idnew[3];
			j=strlen(namefrom); k=strlen(append);
			nameto=(char *)tmalloca(j+k+1,'c');   /* make combined name */
			Loop(i,0,j) nameto[i]=namefrom[i];
			Loop(i,0,k) nameto[i+j]=append[i];
			nameto[j+k]='\0';
			printout("normal"," to %s, sign=",nameto);
			if (af->type=='d') 
			{ v=(double *)createarray(nameto,iallnew*irep,af->type,0);
				vf=(double *)af->pointer;
				Loop(jrep,0,irep)
				{ sign=readint(fpin);
					if (sign!=-1) sign=1;
					printout("normal"," %d",sign);
					Loop(iiold[3],0,i4d[3]) Loop(iiold[L2],0,idold[L2]) 
					Loop(iiold[L3],0,idold[L3])  Loop(iiold[L],0,idold[L])
					{ Loop(i,0,4) ii[i]=iiold[i];
						vv=vf[In4(iiold,idold)+jrep*iallold];
						if (side==0) ii[L]=iiold[L]+idnew[L]-idold[L];  /* copy */
						v[In4(ii,idnew)+jrep*iallnew]=vv;
						ii[L]=idnew[L]-1-ii[L];
						if (iadd==1 && ((side==0 && iiold[L]==0) 
											 || (side!=0 &&  iiold[L]==idold[L]-1))) continue;
						else  v[In4(ii,idnew)+jrep*iallnew]=vv*sign;
					}
				}
				printout("normal","\n");
			}
			else if (af->type=='c') 
			{ c=(char *)createarray(nameto,af->size,af->type,0);
				cf=(char *)af->pointer; 
				Loop(jrep,0,irep)
				Loop(iiold[3],0,i4d[3]) Loop(iiold[L2],0,idold[L2]) 
				Loop(iiold[L3],0,idold[L3])  Loop(iiold[L],0,idold[L])
				{ Loop(i,0,4) ii[i]=iiold[i];
					cc=cf[In4(iiold,idold)+jrep*iallold];
					if (side==0) ii[L]=iiold[L]+idnew[L]-idold[L];  /* copy */
					c[In4(ii,idnew)+jrep*iallnew]=cc;
					ii[L]=idnew[L]-ii[L];
					if (iadd==1 && iiold[L]==idold[L]-1) continue;
					else  c[In4(ii,idnew)+jrep*iallnew]=cc;
				}
			}
		}
		
	}
	tmalloca(-1,'d');
}
