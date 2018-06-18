/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/*  Serie 8 - Tridiagonalmatrizen Inverse & Rayleigh Iteration	 */
/* ------------------------------------------------------------- */
/*	Autoren: 		 		 		 */
/*	Versionsnummer: 					 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "matrix.h"


/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int
main(void){

  int n, i;
  double *x, *al, *ad, *au;
  double norm, lambda, mu;
  double h;

  n = 100;
  h = 1.0 /(1.0*(n+1));

  x = (double*) calloc(n, sizeof(double));
  al = (double*) calloc(n-1, sizeof(double));
  ad = (double*) calloc(n, sizeof(double));
  au = (double*) calloc(n-1, sizeof(double));

  /* ------------------------------------------------------------
   *  Befuellen
   * ------------------------------------------------------------ */   

   for(i = 0; i < n-1; i++){
	 al[i] = -1.0;
     ad[i] = 2.0 / h/ h;
     au[i] = -1.0 /h/ h;
     al[i] = -1.0 /h/ h;		
     x[i] = 2.0 * rand() / RAND_MAX - 1.0;   
   }
     x[n-1] = 2.0 * rand() / RAND_MAX - 1.0;  
   ad[i] = 2.0/ h/ h;
   
 /* ------------------------------------------------------------
   * Inverse Iteration mit Shift
   * ------------------------------------------------------------ */

  printf("Inverse Iteration mit Shift\n");
  lambda = 0.0;
  norm = 0.0;
  mu = 7.0;
  inverse_iteration_withshift_tridiag(n, al, ad, au, x, mu, 10, &lambda, &norm);

  printf("  Eiegenwert %g\n"
	 "  Residuum (Norm) %e\n",
	 lambda, norm);
	 
  /* ------------------------------------------------------------
   * Rayleigh Iteration
   * ------------------------------------------------------------ */

   for(i = 0; i < n; i++){	
     x[i] = 2.0 * rand() / RAND_MAX - 1.0;     
   }


  printf("Rayleigh Iteration\n");
  lambda = 0.0;
  norm = 0.0;
  mu = 7.0;
  rayleigh_iteration_tridiag(n, al, ad, au, x, mu, 10, &lambda, &norm);

  printf("  Eigenwert %g\n"
	 "  Residuum (Norm) %e\n",
	 lambda, norm);


  free(x);
  free(al);
  free(ad);
  free(au);  

  return EXIT_SUCCESS;
}
