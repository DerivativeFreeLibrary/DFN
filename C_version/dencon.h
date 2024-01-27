/* variables, arrays, and routines for dencon */
 /******************** variables ************************/
  int nf_max, num_iter, num_funct;
  int iprint, istop;
  int i_corr, i_corr_fall;

  double alfa_max, alfa_stop;
  double obj;
  double viol, violz, violiniz;
 /********************* arrays **************************/
  double xreal[MAX_VARIABLE+1];
  double lb[MAX_VARIABLE+1];
  double ub[MAX_VARIABLE+1];
  double eps[MAX_CONSTRAINT+1];
  double epsiniz[MAX_CONSTRAINT+1];
  double constr[MAX_CONSTRAINT+1];

  double xfstop[MAX_VARIABLE+1][2*MAX_VARIABLE+1+1];
  double direzioni[MAX_VARIABLE+1][MAX_VARIABLE+1];
  double H[MAX_VARIABLE+1][MAX_VARIABLE+1];
 /*************** storage of x vectors *******************/
struct{
double xreal[MAX_VARIABLE+1];
double obj;
double constr[MAX_CONSTRAINT+1];
double constr_violation;
} typedef Xstore;

Xstore xstore[MAX_XSTORE+1];
int nXstores; /* total number of x vectors stored */
int pointXstore; /* pointer to recently used storage location */
/*eject*/
/************************ routines ************************/

/* denconmain.c */
/* main program */

/* dencon.c */
void dencon(double *x);
void displaysolution(char *label);
double funct(double *x);
void linesearchbox_cont(double *x,double *f_,double *d,
                        double *alfa_,double *alfa_d,
                        double *z1,double *fz1_,
                        double *z2,double *fz2,
                        double *z,double *fz_);
void linesearchbox_dense(double *x,double *f_,double *d,
                        double *alfa_,double *alfa_d_,
                        double *z1,double *fz1_,
                        double *z2,double *fz2,
                        double *z,double *fz_);
void teststop();

/* halton.c */
void halton(int index,double *xx);
void gen_base(double *d1);
void gram_schmidt();


/* problem.c */
void fconstr(double *x);
double fobj(double *x);
void setbounds();
void setdim();
void startp(double *x);

/* xstore.c */
double computeXvalues(double *x);
void retainXstore(double *x);

/************* last record of dencon.h **********/
