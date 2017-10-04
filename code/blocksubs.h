int blockcoef(char m, int *cn, int **ci, double **cc, int *i4dp, int *wherep, 
              char *cltp, int **acnb, int ***acib, double ***accb, int *aneqsb, 
              double **acplus);
void blockrhs(char m, int *i4dp, int *wherep, double *rhs, int *cnb, double **arhsb);
void blocktdma(char m, double *cplus, int *cnb, int **cib, double **ccb, 
               double *rhsb, double **adpb);
void blocktop(char m, int *i4dp, double *dpb, int *cnb, int **cib, double *dp);
