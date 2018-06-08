/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 	Serie 7 - Test Vektor- und Inverse-Iteration	 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 		 				 */
/*	Versionsnummer:						 */
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
  double eps, h, ih;
  double *a, *x;
  time_t t;
  FILE *F; 			/* File */
  
  time(&t);
  srand(t);
  
 /* ----------------------------------------
   Startwerte und Genauigkeit festlegen
   ---------------------------------------- */
  
  printf("Start setting up problem\n");
  n = 20;
  eps = 1e-10;
  
  a = (double *) malloc(sizeof(double) * n*n);
  x = (double *) malloc(sizeof(double) * n);
 
  for(i = 0; i < n; i++){
    x[i] = (rand() % 90)/10.0 + 1.0;
  }	
 
  build_matrix1d(a, n);

 /* ----------------------------------------
   Vektoriteration
   ---------------------------------------- */  
 
  printf("Power iteration\n");   
  power_adaptive(a, n, x, eps);
  
  /* Loesung in Datei schreiben */
  F = fopen("./solution_power.txt", "w");
  h = 1.0/(n+1);
  fprintf(F,"%f %f\n", 0.0, 0.0);
  for(i = 1; i <= n; i++){
    ih = i*h;
    fprintf(F,"%f %f\n", ih, x[i-1]);
  }
  fprintf(F,"%f %f\n", 1.0, 0.0);
  fclose(F);
  
 /* ----------------------------------------
   Problem aufsetzten zweites Mal
   ---------------------------------------- */
   
  for(i = 0; i < n; i++){
    x[i] = (rand() % 90)/10.0 + 1.0;
  }	
  
  build_matrix1d(a, n);
  
/* ----------------------------------------
   Inverse Iteration
   ---------------------------------------- */  
  
  printf("Inverse iteration\n");   
  invit_adaptive(a, n, x, eps);
  
  F = fopen("./solution_invit.txt", "w");
  h = 1.0/(n+1);
  fprintf(F,"%f %f\n",0.0,0.0);
  for(i = 1; i <= n; i++){
    ih = i*h;
    fprintf(F,"%f %f\n", ih, x[i-1]);
  }
  fprintf(F,"%f %f\n", 1.0, 0.0);
  fclose(F);
  
  printf("Finished =)\n");   
  
  /* Speicherfreigabe */
  free(a);
  free(x);
  
  return EXIT_SUCCESS;
}
