/* contains fixdudx */
#include "global.h"

double fixdudx(double (*du)[3])

/* fix du for incompressible continuity, consider 2 and 3d flow */
{
  double cont,d;
  int i,im,ic[3];
  cont=du[0][0]+du[1][1]+du[2][2];
  if (cont==0) return (0);
  d=abs(du[0][0])+abs(du[1][1])+abs(du[2][2]);
  im=0;
  for (i=0;i<3;i++) 
  {
    if (abs(du[i][i])>.001*d) {ic[i]=1; im++;} 
    else ic[i]=0; 
  }
  for (i=0;i<3;i++) du[i][i] -= ic[i]*(cont/im);
  return(cont/d);
} 
