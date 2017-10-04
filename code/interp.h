
int  findex(double a, double *a1d, int id, double *f);      
int find2dinterp(double *x, double (*xp)[2][3], double *ff, int nylimit, int itermax);    
void find1dinterp(double *x, double (*xp)[3], double *ff, int nylimit);   
void interp3d(double *x, int *idim, double *f, double *xout);
int findex3d(double *xx, double **x, int *id, int *ilo, int *ihi, double tol, double *fout);



