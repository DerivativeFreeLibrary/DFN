/* main program */
/***************************************************************
 *   Dencon: Linesearch-based Derivative-free Approach for Nonsmooth
 *   Constrained Optimization
 *
 *   This is an implementation in C of the DFNcon algorithm described
 *       in G. Fasano, G. Liuzzi, S. Lucidi, F. Rinaldi. A Linesearch-
 *       based Derivative-free Approach for Nonsmooth Optimization
 *   Copyright (C) K. Truemper
 *
 *   A Fortran90 implementation by G. Fasano G.Liuzzi, S.Lucidi, 
 *   F. Rinaldi called DFN is available on the website 
 *   http://www.dis.uniroma1.it/~lucidi/DFL
 *
 *   This program is free software: You can redistribute it and/or 
 *   modify it under the terms of the GNU General Public License as 
 *   published by the Free Software Foundation, either version 3 of 
 *   the License, or (at your option) any later version.
 *
 *   For details of the GNU General Public License, see
 *   <http://www.gnu.org/licenses/>.
 *
 *   We do not make any representation or warranty, 
 *   expressed or implied, as to the condition,
 *   merchantability, title, design, operation, or fitness 
 *   of dencon for a particular purpose.
 *
 *   We do not assume responsibility for any errors, 
 *   including mechanics or logic, in the operation 
 *   or use of denconn, and have no liability whatsoever 
 *   for any loss or damages suffered by any user as a result of 
 *   seqpen.
 * 
 *   In particular, in no event shall we be liable
 *   for special, incidental, consequential, or
 *   tort damages, even if we have been advised of 
 *   the possibility of such damages.
 *
 **************************************************************/
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>

# include "basic.h"
# include "dencon.h"

void setnomefun(char* nomefun) ;

/*eject*/
int main () {

  int i;
  double objiniz;
  char nomefun[30];    
  FILE *fid; 
  /* initialize arrays */

  /* define nreal and ncon */
  setdim();

  /* define bounds lb and ub */
  setbounds();

  /* define starting xreal vector */
  startp(xreal);

  /* check that xreal satisfies lower and upper bounds */
  for (i=1; i<=nreal; i++) {
    if ((xreal[i] < lb[i]) || (xreal[i] > ub[i])) {
      printf(
      "\n  denconmain: starting point violates bound constraints\n");
      exit(1);
    }
  }
/*eject*/
  /* initialize counts of function calls and iterations */
  num_funct = 0; 
  num_iter = 0;

  /* evaluate initial solution */

  /* obj for xreal */
  /*obj = fobj(xreal);*/

  /* constr for xreal */
  fconstr(xreal);

  /* choice of starting penalty parameter values */

  /* initialize eps using constr */
  for (i=1; i<=ncon; i++) {
    if (constr[i] < 1.0) {
      eps[i] = 1.e-3;
    } else {
      eps[i] = 1.e-1;
    }
  }

  /* compute viol */
  viol = 0.0;
  for (i=1; i<=ncon; i++) {
    viol += max(0.0,constr[i]);
  }
  funct(xreal);

  /* display initial solution */
  displaysolution("initial");
/*eject*/

  for (i=1; i<=ncon; i++) {
    epsiniz[i] = eps[i];
  }
  objiniz = obj;
  violiniz = viol;

  alfa_stop = 1.e-7;
  nf_max = 20000;
  iprint = 0; /* = 0 if no iteration output */
              /* = 1 or 2 if iteration output */

  /* solve problem */
  printf("\nStart the optimizer");
  dencon(xreal); /* uses search directions */
                 /* plus linesearch along coordinate directions */

  printf("\nTermination condition: ");
  if (istop == 1) {
    printf("convergence\n");
  } else if (istop == 2) {
    printf("max number of evaluations\n");
  } else {
    printf(" dencon: error, unknown termination code = %d",istop);
    exit(1);
  }
/*eject*/
  /* evaluate final solution */

  /* obj, constr, viol for xreal */
  /*obj = fobj(xreal);*/
  /*fconstr(xreal);*/
  funct(xreal);

  viol = 0.0;
  for (i=1; i<=ncon; i++) {
    viol += max(0.0,constr[i]);
  }

  /* display final solution */
  displaysolution("final");
  
  fid = fopen("funct.out","a");
  for (i=num_funct+1; i<=20000; i++){
    if((obj > 1.e+32) || (viol > 1.e-6)){
      fprintf(fid,"            NaN\n");
    }else{
      fprintf(fid,"%15.6e\n",obj);	
    }
  }
  fclose(fid);
  
  /*strcpy(nomefun,"Polak 2 + (n-1) constraints");*/
  setnomefun(nomefun);
  
  printf("\n------------------------");
  /* write LaTeX output line */
  printf("\n %s & %4d & %4d & %6d & %15.6e ",
         nomefun,nreal,ncon,num_funct,objiniz);
  printf("& %15.6e & %15.6e & %15.6e & - \\hline",
         violiniz,obj,viol);
  printf("\n------------------------\n");

  fid = fopen("fort.2","w");
  fprintf(fid,"%s    & %4d & %4d & %6d & %15.6e ",
         nomefun,nreal,ncon,num_funct,objiniz);
  fprintf(fid,"& %15.6e & %15.6e & %15.6e & - \\\\\\hline\n",
         violiniz,obj,viol);
  fclose(fid);

  return 0;

}
/************ last record of denconmain.c **************/
