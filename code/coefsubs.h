void coefcadd(int iw, int *cn, int **ci, double **cc, int nocoefs, int jcoef, int n, int *ipt, double *cadd);
void coefcombine(int iw, int *cn, int **ci, double **cc, int irep, double tol);
void coeforder(int icen, int iw, int *cn, int **ci, double **cc, int jrep);
void coefcenterfirst(int icen, int iw, int *cn, int **ci, double **cc, int jrep);
void coeforderv(int icen, int iw, int *cn, int **ci, double **cc, int jrep, int jorder);

/* in c_coefrhs */
void coefrhs(FILE *fprint, const char *namev, const char *namerhs, int jcoef, char oldnew, double fac);
