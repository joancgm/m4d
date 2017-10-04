void eqnrhs(int ns, int nn, int *cn, int **ci, double **cc, double *ch, 
            double *rhsin, double *rhs);
void eqnrhsbij(int ns, int nn, int *cn, int **ci, double **cc, int npts, int neqs,
               double *beqm, double *dbij, double *beqd, char *nybijeq, double *rhs);
void eqncplus(int ns, int nn, int *cn, double **cc, double *cplus);
void eqnerror(int ns, int nn, double *rhs, double *cplus, double *err, int ierrtot);
int eqnerrorwho(int ns, int nn, double *rhs, double *cplus, char *ny, double *err, 
                int *im, int *ix);
void eqnbndry(int *ipeqexit, int *cn, int **ci, double **cc, double *rhsin, double *rhs, double *ch);
void eqncenter(int ns, int nn, double *cplus, double *rhs, int **ci, double *ch);
void eqncenterbij(int ns, int nn, int neqs, int npts, double *cplus, double *beqd,
                  double *rhs, int **ci, char *nybijeq, double *bij, double *dbij);
