/* tdma structure */
typedef struct tdma_ 
               { int icmax; /* number of variables */
                 int *ic;   /* relative equation indices of variables */
                 double *cc; /* coefficient of ic */
                 double *cm; /* coefficient of ic-1 */
                 double *cp; /* coefficient of ic+1 */
                 struct tdma_ *next; /* address of next Tdma */
               }  Tdma;

Tdma* tdmaalloc(int i);
void tdmafree(Tdma *td);
Tdma *tdmafromcoef(int ns, int ne, int *coefn, int **coefi, double **coefc, double *cplus, int *where);
Tdma *tdmaijkcoef(int L, int ns, int ne, int *i4d, int *coefn, int **coefi, double **coefc, double *cplus, int *where, int *match);
void tdmar(double *cm, double *cc, double *cp, double *bb, double *x, int max);
void tdmasolve(Tdma *td1, int ns, int **coefi, double *cplus, double *bb, double *x);
