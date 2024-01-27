/* routines defining problem elattar */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# include "basic.h"
# include "dencon.h"

/***************************************************************
 *
 * subroutines in this file:
 *       void fconstr(double *x)
 *       void fobj(double *x)
 *       void setbounds()
 *       void setdim()
 *       void startp(double *x)
 **************************************************************/
/*eject*/
/**************************************************************
 *   void fconstr(double *x): defines constr values for given 
 *     problem
 *     all constraints are assumed to be <= 0
 *   uses or modifies: nreal, ncon, constr          
 **************************************************************/
void fconstr(double *x) {

  return;

}
/*eject*/
/**************************************************************
 *   double fobj(double *x): returns obj value for given vector x
 *   uses or modifies: nreal       
 **************************************************************/
double fobj(double *x){

  int i;
  double app[103], f;
  double y[52], t[52];

  num_funct++;

  for(i = 1; i<=51; i++) {
    t[i] = ((double)(i) - 1.0) / 10.0;
    y[i] = exp(t[i])/2.0 - exp(-2.0*t[i]);
    y[i] = y[i] + exp(-3.0*t[i])/2.0;
    y[i] = y[i] + 1.5*exp(-1.5*t[i])*sin(7.0*t[i]);
    y[i] = y[i] +       exp(-2.5*t[i])*sin(5.0*t[i]);
    app[i] = x[1]*exp(-x[2]*t[i])*cos(x[3]*t[i] + x[4]);
    app[i] = app[i] + x[5]*exp(-x[6]*t[i]) - y[i];
  }

  for(i = 52; i<=102; i++) {
    app[i] = -app[i-51];
  }

  f = app[1];
  for(i=2; i<=102; i++) {
    if(app[i] > f) f = app[i];
  }

  return f;
}
/*eject*/
/**************************************************************
 *   void setbounds(): defines lb and ub bounds for given problem
 *   uses or modifies: nreal,lb,ub
 **************************************************************/
void setbounds() {

  int j;

  for (j=1; j<=nreal; j++) {
    lb[j] = -100.0;
    ub[j] =  100.0;
  }

  return;

}
/*eject*/
/**************************************************************
 *   void setdim(): defines number of variables nreal
 *                          number of constraints ncon
 *   uses or modifies: nreal, ncon
 **************************************************************/
void setdim() {

  nreal = 6;
  ncon = 0;

  return;

}
/*eject*/
/**************************************************************
 *   void startp(double *x): defines initial vector x
 *   uses or modifies: nreal
 **************************************************************/
void startp(double *x) {

  x[1] =  2.0;
  x[2] =  2.0;
  x[3] =  7.0;
  x[4] =  0.0;
  x[5] = -2.0;
  x[6] =  1.0;

  return;

}

void setnomefun(char* nomefun) {

  strcpy(nomefun,"Elattar                       ");

  return;
} 
/***************** last record of problem.c ****************/
