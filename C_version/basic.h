/* constants */

/* used only for allocation */
#define MAX_ENTRY 256
#define MAX_CONSTRAINT 300
#define MAX_VARIABLE 300
#define MAX_XSTORE 100

/* general constants */
#define INF 1.0e14
#define TRUE  1
#define FALSE 0

/* basic variables for dencon */

 int nreal; /* number of variables */
 int ncon;  /* number of constraints */

/* max and min definition */
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

/* math library */
double pow(double,double);
double sqrt(double x);
/************** last record of basic.h *************/
