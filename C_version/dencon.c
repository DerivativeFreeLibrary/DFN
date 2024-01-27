/* routines solving problem */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "basic.h"
# include "dencon.h"
/**************************************************************
 *
 * subroutines in this file:
 *       void dencon(double *x)
 *       void displaysolution(char *label)
 *       void funct()
 *       void linesearchbox_cont()
 *       void linesearchbox_dense()
 *       void stop()
 **************************************************************/
/*eject*/
/**************************************************************
 *   void dencon(): main solution subroutine        
 **************************************************************/
void dencon(double *x) {

  int cambio_eps;
  int n, i, j;
  int index_halton;
  int imin, imax;
  int tipo_direzione;

  double dconv[MAX_VARIABLE+1],dnr;
  double z[MAX_VARIABLE+1],
         d[MAX_VARIABLE+1], 
         d_dense[MAX_VARIABLE+1];
  double alfa;
  double alfa_d[MAX_VARIABLE+1],
         alfa_diag[MAX_VARIABLE+1],
         alfa_coord[MAX_VARIABLE+1], 
         alfa_dense[MAX_VARIABLE+1];
  double f, fz;
  double maxeps;
  double fstop[2*MAX_VARIABLE+1+1],
         d1[MAX_VARIABLE+1];
  double fz1, fz2, z1[MAX_VARIABLE+1], z2[MAX_VARIABLE+1];
  double fmin, fmax,soglia;
  double doldalfamin, doldalfamax, rapalfa;

  double dval;
/*eject*/
  /* initialization to suppress compiler warning */
  alfa = 0.0;
  fz1 = 0.0;
  fz2 = 0.0;  
  fz = 0.0;
  obj = 0.0;

  /* define n */
  n = nreal;

  istop = 0;
  for (j=1; j<=2*n+1; j++) {
   fstop[j] = 0.0;
  }
  for (j=1;j<=n; j++) {
    alfa_dense[j] = 0.0;
  }
  /*index_halton = 3581*19;*/
  index_halton = 19;
  /*index_halton = 1000;*/
  soglia = 1.e-3;
/*eject*/
  for (i=1; i<=n; i++) {
    alfa_d[i] = max(1.e-3,min(1.0,fabs(x[i])));
    alfa_diag[i] = alfa_d[i];
    alfa_coord[i] = alfa_d[i];
    alfa_dense[i] += alfa_d[i];
    if (iprint >= 1) {
      printf("\n alfainiz[%d] = %g",i,alfa_d[i]);
      fflush(stdout);
    }
  }

  for (i=1; i<=n; i++) {
    alfa_dense[i] /= (double)n;
  }

  if (n > 1) {
    halton(index_halton,d_dense);
  }

  for (i=1; i<=n; i++) {
    for (j=1; j<=n; j++) {
      direzioni[i][j] = 0.0;
    }
    direzioni[i][i] = 1.0;
    d[i] = 1.0;
  }

  f = computeXvalues(x); /*  also computes obj, constr, viol  */
  for (i=1; i<=2*n+1; i++) {
    fstop[i] = f;
  }
/*eject*/
  i_corr = 1;

  tipo_direzione = 0;

  for (i=1; i<=n; i++) {
    for (j=1; j<=2*n+1; j++) {
      xfstop[i][j] = x[i];
    }
    z[i]=x[i];
  }

  if(iprint >= 2) {
    printf("\n ----------------------------------");
    printf("\n finiz = %g",f);
    for (i=1; i<=n; i++) {
      printf("\n xiniz[%d] = %g",i,x[i]);
    }
    fflush(stdout); 
  }
/*eject*/
  while (1) {

    alfa_max = -INF;
    if (n > 1) {
      for (i=1; i<=n; i++) {
        alfa_max = max(alfa_max,alfa_coord[i]);
        alfa_max = max(alfa_max,alfa_diag[i]);
        alfa_max = max(alfa_max,alfa_dense[i]);
      }
    } else {
      alfa_max = max(alfa_max,alfa_d[1]);
    }

    teststop();
    alfa_max = -INF;
    for (i=1; i<=n; i++) {
      alfa_max = max(alfa_max,alfa_d[i]);
    }   

    if (istop >= 1) {
      return;
    }

    if (i_corr == 1) {
      for (i=1; i<=n; i++) {
        dconv[i] = 0.0;
      }
      for (i=1; i<=n; i++) {
        for (j=1; j<=n; j++) {
          dconv[j] -= direzioni[j][i];
        }
      }
    }
/*eject*/
    if (num_iter%20 == 0) {
      printf("\n num_funct = %3d num_iter = %d obj = %g alfamax = %g",
             num_funct,num_iter,obj,alfa_max);
    }

    if (iprint >= 1) {
      printf("\n ----------------------------------------------");
      printf("\n num_iter = %d num_funct = %d f = %g alfa_max = %g",
             num_iter,num_funct,f,alfa_max);
      fflush(stdout);
    }
    if (iprint >= 2) {
      for (i=1; i<=n; i++) {
        printf("\n x[%d] = %g",i,x[i]);
      }
      fflush(stdout);
    }

    for (i=1; i<=n; i++) {
      d[i] = direzioni[i][i_corr];
    }
    if (tipo_direzione == 0) {
      linesearchbox_cont(x,&f,d,&alfa,alfa_d,
                         z1,&fz1,z2,&fz2,
                         z,&fz);
      if(fabs(alfa) >= 1.e-12) {
        x[i_corr] += alfa*d[i_corr];
      }
    } else {
      linesearchbox_dense(x,&f,d,&alfa,&alfa_d[i_corr],
                          z1,&fz1,z2,&fz2,
                          z,&fz);
      if(fabs(alfa) >= 1.e-12) {
        for (i=1; i<=n; i++) {
          x[i] = min(ub[i], x[i]+alfa*d[i]);
          x[i] = max(lb[i],x[i]);
        }
      }
    }

    for (i=1; i<=n; i++) {
      direzioni[i][i_corr] = d[i];
    }
/*eject*/
    if (fabs(alfa) >= 1.e-12) {
      f = fz;
      viol = violz;
      num_iter++;
    } else {

      if(i_corr_fall == 0) { 
        fstop[i_corr] = fz1;
        fstop[2*i_corr] = fz2;
             
        for (j=1; j<=n; j++) {
          xfstop[j][i_corr] = z1[j];
          xfstop[j][2*i_corr] = z2[j];
        }            

        num_iter++;
      }
    }

    for (i=1; i<=n; i++) {
      z[i] = x[i];
    }
/*eject*/
    if (i_corr < n) {
      i_corr++;
    } else {
      dval = -INF;
      for (i=1; i<=n; i++) {
        dval = max(dval,alfa_d[i]);
      }
      if ((dval <= soglia) && (n > 1)) {

        if (tipo_direzione == 0) {

          fmin = fstop[1]; /* = f */
          fmax = fstop[1]; /* = f */
          imin = 1;
          imax = 1;
          doldalfamin = alfa_d[1];
          doldalfamax = alfa_d[1];
          for (i=2; i<=n; i++) {
            if (alfa_d[i] < doldalfamin) {
              doldalfamin = alfa_d[i];
            } 
            if (alfa_d[i] > doldalfamax) {
              doldalfamax = alfa_d[i];
            } 
          }
/*eject*/  
          rapalfa = 3.0;

          if (doldalfamax/doldalfamin > rapalfa) {
          
            for (i=1; i<=n; i++) {
              d1[i] = dconv[i];
            }
            dnr = sqrt((double)n);
          } else {              
            for (i=2; i<=2*n; i++) {
              if (fstop[i] < fmin) {
                 fmin = fstop[i];
                 imin = i;
              }
              if (fstop[i] > fmax) {
                fmax = fstop[i];
                imax = i;
              }
            }

            dnr = 0.0;
            for (i=1; i<=n; i++) {             
              d1[i] = xfstop[i][imin] - xfstop[i][imax];
              dnr += d1[i]*d1[i];
            }
            dnr = sqrt(dnr);

            if (dnr <= 1.e-24) {
              for (i=1; i<=n; i++) {
                 d1[i] = dconv[i];
              }
              dnr = sqrt((double)n);
            }
          }
/*eject*/
          /* define H[][] matrix from d1 vector */
          gen_base(d1);
          /* perform Gram-Schmidt orthogonalization on H[][] */
          gram_schmidt();
          /* transfer H[][] matrix into direzione[][] matrix */
          for (i=1; i<=n; i++) {
            for (j=1; j<=n; j++) {
              direzioni[i][j] = H[i][j];
            }
          }

          tipo_direzione = 1;
          for (i=1; i<=n; i++) {                  
            alfa_coord[i] = alfa_d[i];
          }        
          if (doldalfamax/doldalfamin > rapalfa) {
            for (i=1; i<=n; i++) {
              alfa_d[i] = 1.0*alfa_diag[i];
            }
          } else {
            dnr = 0.0;
            for (i=1; i<=n; i++) { 
              dnr += alfa_d[i];
            }         
            dnr /= (double)n;
            for (i=1; i<=n; i++) {
              alfa_d[i] = 1.0*dnr;
            }
          }
             
          if (iprint >= 1) {
            printf("\n fine dir. coordinate");
            fflush(stdout);
          }
/*eject*/           
        } else if (tipo_direzione == 1) {
          /* define H[][] matrix from d_dense vector */
          gen_base(d_dense);
          /* perform Gram-Schmidt orthogonalization on H[][] */
          gram_schmidt();
          /* transfer H[][] matrix into direzione[][] matrix */
          for (i=1; i<=n; i++) {
            for (j=1; j<=n; j++) {
              direzioni[i][j] = H[i][j];
            }
          }
          /*index_halton +=3581;*/
          index_halton++;
          /*index_halton +=2*n;*/
          halton(index_halton,d_dense);

          tipo_direzione = 2;
          for (i=1; i<=n; i++) {
            alfa_diag[i] = alfa_d[i];
          }
              
          dnr = 0.0;
          for (i=1; i<=n; i++) {
            dnr += alfa_d[i];
          }          
          dnr /= (double)n;

          for (i=1; i<=n; i++) {
            alfa_d[i] =1.e+1*dnr;
          } 

          if (iprint >= 1) {
            printf("\n fine dir. n+1");
            fflush(stdout);
          }
/*eject*/
        } else if (tipo_direzione == 2) {

          for (i=1; i<=n; i++) {
            for (j=1; j<=n; j++) {
              direzioni[i][j] = 0.0;
            }
            direzioni[i][i] = 1.0;
          }

          tipo_direzione = 0;

          for (i=1; i<=n; i++) {
            alfa_dense[i] = alfa_d[i];
            alfa_d[i] = 1.0*alfa_coord[i];
          }         

          if (iprint >= 1) {
            printf("\n fine dir. densa");
            fflush(stdout);
          }

        } else {
          printf("\n dencon: tipo_direzione = %d != 0, 1, 2",
                 tipo_direzione);
          exit(1);
        } /* end if tipo_direzione == 0, else if ... */
      
        i_corr = 1; 

      } /* end if dval < soglia ... */
   
      i_corr = 1;

    } /* end if i_corr < n, else */
/*eject*/   
    f = computeXvalues(x); /* also computes obj, constr, viol */

    if ((viol > 0.0) && (ncon >= 1)) {      
      cambio_eps = FALSE;
      maxeps = -INF;
      for (i=1; i<=ncon; i++) {
        maxeps = max(maxeps,eps[i]);
      }
      dval = -INF;
      for (i=1; i<=n; i++) {
        dval = max(dval,alfa_d[i]);
      }

      if (iprint >= 1) {
        printf("\n**************************************");
        printf("\n**************************************");
      }
      for (i=1; i<=ncon; i++)  {
        if (eps[i]*constr[i] > 1.0*dval) {
          eps[i] *= 1.e-2;
          if (iprint >= 1) {
            printf(
            "\n*********** aggiorno eps[%d] = %g *************",
            i, eps[i]);
            fflush(stdout);
          }
          cambio_eps = TRUE;
        }
      }
      if (iprint >= 1) {
        printf("\n**************************************");
        printf("\n**************************************");
      }
/*eject*/
      if (cambio_eps == TRUE) {

        /* recompute f using new eps[] */
        f = obj;
        for (i=1; i<=ncon; i++) {
          if (eps[i] <= 0.0) {
            printf("\n dencon: error, eps[%d] = %g <= 0.0",i,eps[i]);
            exit(1);
          }
          f += max(0.0,constr[i])/eps[i];
        }

        /* recompute alfa_d[] */
        for (i=1; i<=n; i++) {
          alfa_d[i] = max(1.e-3,min(1.0,fabs(x[i])));
          if (iprint >= 1) {
            printf("\n alfainiz[%d] = %g",i, alfa_d[i]);
            fflush(stdout);
          }
        } 
      }
    }

  } /* end while (1) */

  return;

}
/*eject*/
/**************************************************************
 *   void displaysolution(char *label): display
 *   obj, xreal, constr, eps, viol
 *   under the heading given by label string
 **************************************************************/
void displaysolution(char *label) {

  int i;

  printf("\n------------------------");
  printf("\n---- %s values ----",label);
  printf("\n------------------------");
  /* obj */
  printf("\nobj = %g\n",obj);
  /* xreal */
  for (i=1; i<=nreal; i++) {
    printf("\nxreal[%d] = %g",i,xreal[i]);
  }
  printf("\n");
  /* constr and eps */
  for (i=1; i<=ncon; i++) {
    printf("\nconstr[%d] = %g\teps[%d] = %g",
           i,constr[i],i,eps[i]);
  }
  printf("\n");
  /* viol */
  printf("\nmax violation value = %g",viol);
  if (viol == 0.0) {
    printf("\nsolution is feasible");
  } else {
    printf("\nsolution is infeasible");
  }
/* num_iter and num_funct */
  if (strcmp(label,"final") == 0) {
    printf("\nnumber of iterations = %d",num_iter);
    printf("\nnumber of obj/constr evaluations = %d",num_funct);
  } 
  printf("\n------------------------");
  printf("\n");
  fflush(stdout);

  return;

}
/*eject*/
/**************************************************************
 *   double funct(double *x): evaluate obj function and 
 *          constraints and return combined measure fmax
 *          also computes obj value obj and violation value viol
 **************************************************************/
double funct(double *x) {

  int i;
  double fmax;
  FILE* fid;

  /* compute obj; initialize fmax and viol */
  obj = fobj(x);
  fmax = obj;
  viol = 0.0;
  /* compute constr for x */
  fconstr(x);
  /* update fmax and viol using constr and eps */
  for (i=1; i<=ncon; i++) {
    if (eps[i] <= 0.0) {
      printf("\n funct: error, eps[%d] = %g <= 0.0",i,eps[i]);
      exit(1);
    }
    fmax += max(0.0,constr[i])/eps[i];
    viol += max(0.0,constr[i]);
  }

  fid = fopen("funct.out","a");
  if((obj > 1.e+32) || (viol > 1.e-6)){
    fprintf(fid,"            NaN\n");
  }else{
    fprintf(fid,"%15.6e\n",obj);	
  }
  fclose(fid);

  return fmax;

}
/*eject*/
/**************************************************************
 *   void linesearchbox_cont(double *x,double *f_,double *d,
 *                      double *alfa_,double alfa_d,
 *                      double *z1,double *fz1_,double *z2,double *fz2,
 *                      double *z,double *fz_)
 *     linearsearch process 'cont'         
 **************************************************************/
void linesearchbox_cont(double *x,double *f_,double *d,
                       double *alfa_,double *alfa_d,
                       double *z1,double *fz1_,double *z2,double *fz2_,
                       double *z,double *fz_) {
  /* transfer variables */
  double f, alfa, fz1, fz2, fz;

  /* local variables */
  int i,j, n;
  int ifront,ielle;
  double alfaex, gamma;
  double delta,delta1,fpar,fzdelta,violzdelta;


  /* define n */
  n = nreal;

  /* transfer in */
  f = *f_;
  alfa = *alfa_;
  fz1 = *fz1_;
  fz2 = *fz2_;
  fz = *fz_;
/*eject*/
  gamma = 1.e-6;   
  delta = 0.5;
  delta1 = 0.5;

  i_corr_fall = 0;
  ifront = 0;
  j = i_corr;

  if (iprint >= 1) {
    printf("\n variabile continua  j = %d    d[j] = %g alfa = %g",
           j, d[j], alfa_d[j]);
    fflush(stdout);
  }


  if (fabs(alfa_d[j]) <= 1.e-3*min(1.0,alfa_max)) {
    alfa = 0.0;
    if (iprint >= 1) {
       printf("\n alfa piccolo");
       printf("\n alfa_d[j] = %g    alfamax = %g",
              alfa_d[j], alfa_max);
       fflush(stdout);
    }
    goto transfer_out;
  }
/*eject*/     
  for (ielle=1; ielle<=2; ielle++) {

    if (d[j] > 0.0) {

      if ((alfa_d[j]-(ub[j]-x[j])) < -1.e-6) {                 
          alfa = max(1.e-24,alfa_d[j]);
      } else {
        alfa = ub[j]-x[j];
        ifront = 1;
        if (iprint >= 1) {
          printf("\n punto espan. sulla front. *");
          fflush(stdout);
        }
      }

    } else {

      if ((alfa_d[j]-(x[j]-lb[j])) < -1.e-6) {
        alfa = max(1.e-24,alfa_d[j]);
      } else {
        alfa = x[j]-lb[j];
        ifront = 1;
        if (iprint >= 1) {
          printf("\n punto espan. sulla front. *");
          fflush(stdout);
        }
      }

    }
/*eject*/
    if(fabs(alfa) <= 1.e-3*min(1.0,alfa_max)) {
  
      d[j] = -d[j];
      i_corr_fall++;
      ifront=0;

      if(iprint >= 1) {
        printf("\n direzione opposta per alfa piccolo");
        printf("\n j = %d    d[j] = %g", j, d[j]);
        printf("\n alfa = %g    alfamax = %g", alfa, alfa_max);
        fflush(stdout);
      }
      alfa = 0.0;
      continue;

    }
/*eject*/
    alfaex = alfa;

    z[j] = x[j] + alfa*d[j];   
    fz = computeXvalues(z); /*  also computes obj, constr, viol  */

    if (ielle == 1) {
      for (i=1; i<=n; i++) {
        z1[i] = z[i];
      }
      fz1 = fz;
    } else {
      for (i=1; i<=n; i++) {
        z2[i] = z[i];
      }
      fz2 = fz;
    }

    violz=viol;
    if (iprint >= 1) {
      printf("\n fz = %g   alfa = %g", fz, alfa);
      fflush(stdout);
    }
    if (iprint >= 2) {
      for (i=1; i<=n; i++) {
        printf("\n z[%d] = %g",i,z[i]);
        fflush(stdout);
      }
    }

    fpar = f - gamma*alfa*alfa;
    if (fz < fpar) {
/*eject*/
      while (1) { /* while #1 */

        if((ifront == 1)) {
          if(iprint >= 1) {
            printf(
              "\n accetta punto sulla frontiera fz = %g   alfa = %g",
              fz, alfa);
            fflush(stdout);
          }
          alfa_d[j] = delta * alfa;
          goto transfer_out;
        }

        if (d[j] > 0.0) {
              
           if (alfa/delta1-(ub[j]-x[j]) < -1.e-6) {
             alfaex = alfa/delta1;
           } else {
             alfaex = ub[j] - x[j];
             ifront = 1;
             if(iprint >= 1) {
              printf("\n punto espan. sulla front.");
              fflush(stdout);
             }
           }
/*eject*/
         } else {

           if(alfa/delta1-(x[j]-lb[j]) < -1.e-6) {
             alfaex = alfa/delta1;
           } else {
             alfaex = x[j]-lb[j];
             ifront = 1;
             if(iprint >= 1) {
              printf("\n punto espan. sulla front.");
              fflush(stdout);
             }
           }

         }
             
         z[j] = x[j] + alfaex*d[j]; 
                   
         fzdelta = computeXvalues(z); 
            /* also computes obj, constr, viol  */
         violzdelta = viol;

         if(iprint >= 1) {
            printf("\n  fzex = %g  alfaex = %g",fzdelta, alfaex);
            fflush(stdout);  
         }
         if(iprint >= 2) {
           for (i=1; i<=n; i++) {
             printf("\n z[%d] = %g",i,z[i]);
             fflush(stdout);
           }
         }
/*eject*/
         fpar = f - gamma*alfaex*alfaex;

         if (fzdelta < fpar) {

           fz = fzdelta;
           violz = violzdelta;
           alfa = alfaex;

         } else {               

           alfa_d[j] = delta*alfa;
           if(iprint >= 1) {
             printf("\n accetta punto fz = %g   alfa = %g", fz, alfa);
             fflush(stdout);
           }
           goto transfer_out;
         }

       } /* end while #1 */
/*eject*/
     } else {      

       d[j] = -d[j];
       ifront = 0;

       if(iprint >= 1) {
         printf("\n direzione opposta");
         printf("\n j = %d    d[j] = %g",j, d[j]);
         fflush(stdout);
       }

     }            
  } /* end for ielle=1 */

  if (i_corr_fall != 2) {
    alfa_d[j] = delta * alfa_d[j];
  }

  alfa = 0.0;

  if (iprint >= 1) {
    printf("\n fallimento direzione");
    fflush(stdout);
  }
/*eject*/
transfer_out:;
  *f_ = f;
  *alfa_ = alfa;
  *fz1_ = fz1;
  *fz2_ = fz2;
  *fz_ = fz;

  return;

}
/*eject*/
/**************************************************************
 *   void linesearchbox_dense(double *x,double *f_,double *d,
 *                      double *alfa_,double *alfa_d_,
 *                      double *z1,double *fz1_,double *z2,double *fz2,
 *                      double *z,double *fz_)
 *     linearsearch process 'dense'         
 **************************************************************/
void linesearchbox_dense(double *x,double *f_,double *d,
                       double *alfa_,double *alfa_d_,
                       double *z1,double *fz1_,double *z2,double *fz2_,
                       double *z,double *fz_) {
  /* caution: alfa_d here is a scalar, not a vector */
  /* transfer variables */
  double f, alfa, alfa_d, fz1, fz2, fz;

  /* local variables */
  int i, n;
  int ielle;
  double alfaex, gamma;
  double delta,delta1,fpar,fzdelta,violzdelta;


  /* define n */
  n = nreal;

  /* transfer in */
  f = *f_;
  alfa = *alfa_;
  fz1 = *fz1_;
  fz2 = *fz2_;
  fz = *fz_;
  alfa_d = *alfa_d_;
/*eject*/
  gamma = 1.e-6;    

  delta = 0.5;
  delta1 = 0.5;

  i_corr_fall = 0;
 
  if (iprint >= 1) {
    printf("\n direzione halton, alfa = %g",alfa_d);
    fflush(stdout);
  }
/*eject*/
  for (ielle=1; ielle<=2; ielle++) {

    alfa = alfa_d;
    alfaex = alfa;
    for (i=1; i<=n; i++) {
      z[i] = x[i] + alfa*d[i];
    }
    for (i=1; i<=n; i++) {
      z[i] = max(lb[i],min(ub[i],z[i]));
    }
   
    fz = computeXvalues(z); /*  also computes obj, constr, viol  */

    if (ielle == 1) {
      for (i=1; i<=n; i++) {
        z1[i]  = z[i];
      }
        fz1 = fz;
    } else {
      for (i=1; i<=n; i++) {
        z2[i]  = z[i];
      } 
      fz2 = fz;
    }

    violz = viol;
/*eject*/
    if (iprint >= 1) {
      printf("\n fz = %g   alfa = %g",fz, alfa);
      fflush(stdout);
    }
    if (iprint >= 2) {
      for (i=1;i<=n; i++) {
        printf("\n z[%d] = %g",i, z[i]);
        fflush(stdout);
      }
    }

    fpar = f - gamma*alfa*alfa;

    if (fz < fpar) {
/*eject*/
      while (1) {

        alfaex = alfa/delta1;

        for (i=1; i<=n; i++) {          
          z[i] = x[i] + alfaex*d[i];          
          z[i] = max(lb[i],min(ub[i],z[i]));
        }    
         
        fzdelta = computeXvalues(z); 
           /* also computes obj, constr, viol  */
        violzdelta = viol;

        if (iprint >= 1) {
          printf("\n fzex = %g  alfaex = %g",fzdelta,alfaex);
          fflush(stdout);  
        }
        if(iprint >= 2) {
          for (i=1; i<=n; i++) {
            printf("\n z[%d] = %g",i, z[i]);
            fflush(stdout);
          } 
        }

        fpar = f - gamma*alfaex*alfaex;

        if (fzdelta < fpar) {

          fz = fzdelta;
          violz = violzdelta;
          alfa = alfaex;

        } else {               
          alfa_d = alfa;
          if (iprint >= 1) {
            printf("\n denso: accetta punto fz = %g   alfa = %g",
                   fz, alfa);
            fflush(stdout);
          }
          goto transfer_out;
        }

      }  /* end while (1) */
/*eject*/
    } else {      

      for (i=1; i<=n; i++) {
        d[i] = -d[i];
      }

      if (iprint >= 1) {
        printf("\n denso:  direzione opposta");
      }

    }  
        
  } /* end for ielle=1 */     

  alfa_d = delta*alfa_d;
  alfa=0.0;

  if (iprint >= 1) {
    printf("\n denso: fallimento direzione");
    fflush(stdout);
  }
/*eject*/
transfer_out:;
  /* transfer out */
  *f_ = f;
  *alfa_ = alfa;
  *fz1_ = fz1;
  *fz2_ = fz2;
  *fz_ = fz;
  *alfa_d_ = alfa_d;

  return;

}
/*eject*/
/**************************************************************
 *   void teststop(): evaluate stopping criteria          
 **************************************************************/
void teststop() {

  istop = 0;

  if (alfa_max <= alfa_stop) {
    istop = 1;
  }
      
  if (num_funct > nf_max) {
    istop = 2;
  }

  return;

}
/********* last record of dencon.c **********/
