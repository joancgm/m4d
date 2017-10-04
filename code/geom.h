/* types for geom8 and geomc functions */
void geom8init(void);
void geom8volfmid(double *vol, double (*fmid)[3], int *ip);
void geom8fx27(double (*f)[3], double (*xg)[3], int *ip);
void geom8areamid(double *area, double *xmid, int *id, int L);
void geom8gradvol(double (*gradv)[8][3], int *ip);
void geom8gradvxf(double (*gradv)[8][3], double (*x)[3], double (*f)[3]);
void geom8vol(double *vol, int *ip);
void geomcinit(void);
void geomcavector(double (*area)[3], int *ic, int L);
void geomcvave(double *cg, int *ip);
double geomcgradvol(int *ip, double (*grad)[3]);
