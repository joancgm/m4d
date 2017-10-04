/* header for plot routines */

/* plotsubs  */
void p_outjgm(char*filename, int*image,int*ipix);
void p_setpixel(int i,int j,int pt_color,int *image,int *ipix);
void p_stline(int i1,int i2,int j1,int j2,int line_width,int line_color,
              int *image,int *ipix);
void p_stlinex(double x1, double x2, double y1,double y2,int line_width,int line_color,
               int *image,int *ipix);
void p_stlimited(double xp1,double yp1,double xp2,double yp2,int linew,int ic, 
                 double xb1, double yb1,double xb2,double yb2, int *image);
void p_fillx2d(double x1, double x2, double x3, double y12, double y3,
               int ic, int *image,int *ipix);
void p_filltri2d(double x1, double x2, double x3, double y1, double y2, double y3,
                 int ic, int *image,int *ipix);
double p_xlin(double x1, double y1, double x2, double y2, double y);
void p_painttri(double *xy[3],double p[3],char nyout, 
              int* ic, double *pc, int ncolor, int *image, int ipix[2]);
double p_xlin(double x1, double y1, double x2, double y2, double y);

/* m4dplot specific commands listed in commandlist.h */

/* outgif */
void p_outgif(char*filename, int*image, int*ipix);
int p_newline(FILE *fv);
void p_set2char(int i,char c2[2]);
void p_gifheader(FILE *fg,int ixtot,int iytot,int numpix,char *rgb);
void p_writebits(int iv, int nbits, FILE *fp, int *size, int *nextbit, int outbuf[200]);
void p_lzw(FILE *fp,int numpix, int imax, int* ic);

/* p_symbol */
FILE *p_openfile(char *name);
int p_readfonts(int debug);
void p_symbolx(double x, double y,char *string,int ifont,int leftright,int updown,
               int aorp,int color, int* imagetot);
void p_symbol(int ix,int iy,char *string,int ifont,int leftright,int updown,
              int aorp,int color, int* imagetot); 

/* other */
void p_glines(int *gcolor,int linewidth,
         int*i4d, double *x[3],char *clt,
         int irange[3][2],double *gxx[2],double scale,double offset[2],
         int *image,int ipix[2]);
void p_vectors(double *u[3], double *vparms,
               int*i4d, double *x[3],char *clt,
               int irange[3][2],double *gxx[2],double scale,double offset[2],
               int *image, int ipix[2]);
void p_fillcontour(double *p, int nystep, int iadd, char ** fnames, double *fparms,
                   int*i4d, double *x[3],char *clt,
                   int irange[3][2],double *gxx[2],double scale,double offset[2],
                   int *image, int ipix[2]);

double p_nextpt2d(double dxy[2], double *x[2],int id[2], 
                  int ilowin[2], double fin[2], int ilowout[2], double fout[2]);
