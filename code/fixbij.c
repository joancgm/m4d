/* contains fixbij fixdbij */
#include "global.h"
/*     fix bij at a point for realizability  */

double fixbij(double *bij, double gtol)
{  
   int i,j,k,itmax=11;
   double b[3][3],bb,bbb,g,gmin,gmax,f,fmin,fmax;
	
   f=1;
   b[0][0]=bij[0];
   b[1][1]=bij[1];
   b[2][2]=bij[2];
   b[0][1]=bij[3]; b[1][0]=b[0][1];
   b[0][2]=bij[4]; b[2][0]=b[0][2];
   b[1][2]=bij[5]; b[2][1]=b[1][2];
   bb=0.; bbb=0.;
   for (i=0;i<3;i++) for (j=0;j<3;j++)
	{ 
     bb += b[i][j] * b[i][j];
		for (k=0;k<3;k++) bbb += b[i][j] * b[j][k] * b[k][i];
	}
   g=1.-4.5*bb+9.*bbb;
	if (bb<=2./3. && g>=0) return f;
	
	fmin=0; fmax=1.; gmin=1.; gmax=g;
	if (bb>2./3.) 
	{ 
     fmax=sqrt(2./3./bb);
		gmax=1-fmax*fmax*(4.5*bb+9*bbb*fmax);
	}
	for (j=0;j<itmax;j++)  /* alternate interp and averaging look for g=0 */
	{  
		if (j%2==0) f=fmin+(fmax-fmin)*(0-gmin)/(gmax-gmin);
		else f=.5*(fmin+fmax);
		g=1-f*f*(4.5*bb-9*bbb*f);
		if (g>0) {fmin=f; gmin=g;}
		else {fmax=f; gmax=g;}
		if (g>=0 && g< gtol) break;
	}
	if (g<0) {f=fmin; g=gmin;}
	for (i=0;i<6;i++) bij[i]*=f;
		return f;
}
/* ------------------------------- */
double fixdbij(double *bij, double *dbij, double gtol)
{  
   int i,j,k,itmax=11;
   double b[3][3],d[3][3];
	double bb=0,bd=0,dd=0,bbb=0,bbd=0,bdd=0,ddd=0;
	double b2,b3,g,f,fn,fmin,fmax;
	
   f=1;
   b[0][0]=bij[0];
   b[1][1]=bij[1];
   b[2][2]=bij[2];
   b[0][1]=bij[3]; b[1][0]=b[0][1];
   b[0][2]=bij[4]; b[2][0]=b[0][2];
   b[1][2]=bij[5]; b[2][1]=b[1][2];
   d[0][0]=dbij[0];
   d[1][1]=dbij[1];
   d[2][2]=dbij[2];
   d[0][1]=dbij[3]; d[1][0]=d[0][1];
   d[0][2]=dbij[4]; d[2][0]=d[0][2];
   d[1][2]=dbij[5]; d[2][1]=d[1][2];
	
   for (i=0;i<3;i++) for (j=0;j<3;j++)
	{ 
     bb+= b[i][j]*b[i][j];
		bd+= b[i][j]*d[i][j];
		dd+= d[i][j]*d[i][j];
		for (k=0;k<3;k++)
		{ 
        bbb += b[i][j] * b[j][k] * b[k][i];
			bbd += b[i][j] * b[j][k] * d[k][i];
			bdd += b[i][j] * d[j][k] * d[k][i];
			ddd += d[i][j] * d[j][k] * d[k][i];
		}
	}
	b2=bb+dd+2*bd;
	b3=bbb+ddd+3*(bbd+bdd);
   g=1.-4.5*b2+9.*b3;
	if (b2<=2./3. && g>=0) return f;  /* ok */
	
	fmin=0; fmax=1;
	g=1-4.5*bb+9*bbb;
	if (bb>2./3. || g<0) /* bij not realizable, correct that and return dbij=0 */
	{
		f=fixbij(bij,gtol); 
		for (i=0;i<6;i++) dbij[i]=0;
		f=-1+f; return f;
	}
	fn=-bd/dd;  /* try first min in b2 if in range 0 to 1*/
	if (fn<=0 || fn>=1) fn=.5*(fmin+fmax);
		
	for (j=0;j<itmax;j++)
		{ 
        f=fn;
			b2=bb+f*(2*bd+f*dd);
			b3=bbb+f*(3*bbd+f*(3*bdd+f*ddd));
			g=1-4.5*b2+9*b3;
			if (b2>2./3. || g<0) {fmax=f;}
			else 
			{ 
           fmin=f;
				if (g<gtol) break;
			}
			fn=.5*(fmin+fmax);
		}
	if (g<0 || b2>2./3.) f=fmin;
		for (i=0;i<6;i++) dbij[i]*=f; 
			return f;
}
