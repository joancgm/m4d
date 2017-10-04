/* this file is not included, but where needed individual macros copied*/

#define Cross(a,b,c) *(c)= *((a)+1) * *((b)+2) - *((a)+2) * *((b)+1); *((c)+1)= *((a)+2) * *(b) - *(a) * *((b)+2);  *((c)+2)= *(a) * *((b)+1) - *((a)+1) * *(b) 
#define Dot(a,b) ( *(a) * *(b) + *((a)+1) * *((b)+1) + *((a)+2) * *((b)+2) )
#define Det(a,b,c) (*(c) *(*((a)+1) * *((b)+2) - *((a)+2) * *((b)+1)) +  *((c)+1) *( *((a)+2) * *(b) - *(a) * *((b)+2)) +   *((c)+2) *( *(a) * *((b)+1) - *((a)+1) * *(b) ))
#define In4(ii,idim) (ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+idim[2]* ii[3])))
#define Loop(n,a,b) for (n=a;n<b;n++)
#define Loop3(ii,a,idim) Loop(ii[0],a,idim[0]) Loop(ii[1],a,idim[1]) Loop(ii[2],a,idim[2])
#define Prod4(ii) (ii[0]*ii[1]*ii[2]*ii[3])
#define Value3d(x,idim,ii,ff) ((1.-ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*ii[2])] +  (ff[0])*(1.-ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*ii[2])] +  (1.-ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*ii[2])] +  (ff[0])*(ff[1])*(1.-ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*ii[2])] + (1.-ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (ff[0])*(1.-ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+idim[1]*(ii[2]+1))] +  (1.-ff[0])*(ff[1])*(ff[2])*x[ii[0]+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))] +  (ff[0])*(ff[1])*(ff[2])*x[ii[0]+1+idim[0]*(ii[1]+1+idim[1]*(ii[2]+1))])
#define Sqr(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define Vabs(a) (sqrt((double)( *(a) * *(a) + *((a)+1) * *((a)+1) + *((a)+2) * *((a)+2) )))
#define Limit(a,b,c) (max(min( (a) ,(c) ),(b) ) )