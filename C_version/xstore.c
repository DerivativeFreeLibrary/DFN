/* storage and retrieval of solutions and associated values */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "basic.h"
# include "dencon.h"
/**************************************************************
 *
 * subroutines in this file:
 *       void computeXvalues(double *x)
 *       double retrieveXstore(double *x)
 **************************************************************/
/*eject*/
/**************************************************************
 *   double computeXvalues(double *x): 
 *          retrieve the following values from xstore array 
 *          for input vector x, or compute directly with funct():
 *              obj, constr[], viol
 *          return fmax = obj + penalty sum using current eps[]
 **************************************************************/
double computeXvalues(double *x) {

  int flag, i, k;
  double fmax;

  for (k=1; k<=nXstores; k++) {
    flag = TRUE;
    for (i=1; i<=nreal; i++) {
      if (x[i] != xstore[k].xreal[i]) {
        flag = FALSE;
        break;
      } 
    }    
    if (flag == TRUE) {
      /* have found xstore vector matching input vector x */
      obj = xstore[k].obj;
      for (i=1; i<=ncon; i++) {
        constr[i] = xstore[k].constr[i];
      }
      viol = xstore[k].constr_violation;
      /* compute fmax using current eps[] values */
      fmax = obj;
      for (i=1; i<=ncon; i++) {
        fmax += max(0.0,constr[i])/eps[i];
      }
      return fmax;
    }
  }

  /* x is not in Xstore */
  /*  compute values with funct() and store in xstore[] */
  fmax = funct(x);
  retainXstore(x);
  return fmax;

}
/*eject*/
/**************************************************************
 *   retainXstore(double *x):
 *        retain input vector x, obj, constr, viol in xstore[]
 **************************************************************/
void retainXstore(double *x) {

  int i;

  if (nXstores < MAX_XSTORE) {
    nXstores++;
    pointXstore = nXstores;
  } else {
    /* storage is full; advance pointer */
    pointXstore++;
    if (pointXstore > MAX_XSTORE) {
    /* use pointer in circular fashion */
      pointXstore = 1;
    }
  }
  for (i=1; i<=nreal; i++) {
    xstore[pointXstore].xreal[i] = x[i];
  }
  xstore[pointXstore].obj = obj;
  for (i=1; i<=ncon; i++) {
    xstore[pointXstore].constr[i] = constr[i];
  }
  xstore[pointXstore].constr_violation = viol;

  return;

}
/********** last record of xstore.c **************/
