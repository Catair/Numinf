/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/*	Serie 5 - Loesen mit QR-Zerlegung			 */
/* ------------------------------------------------------------- */
/*	Autoren: 				 		 */
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

int main(void){

  /* ---------------------------------------------- */
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */
  /*     T       O     O     D     D     O     O    */
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */
  /* ---------------------------------------------- */

  return EXIT_SUCCESS;
}

double
get_entry(double* a, int ld, int row, int col){

    return a[row+ld*col];

}

void
set_entry(double* a, int ld, int row, int col, double value){

    a[row+ld*col] = value;

}

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void
backward_subst(int ld, double* r, double* b, double* x){

    int i,j;

    for(j=(ld-1);j>=0;j--) {
        x[j] = b[j] / r[(ld+1)*j];
        for(i=0;i<j;i++) {
            b[i]= b[i]-r[i+j*ld]*x[j];
        }
    }
}

void
print_matrix(double* a, int rows, int cols){
    int i, j;

    for(j=0;j<rows;j++) {
        printf(" | ");

        for(i=0;i<cols;i++) {

            printf(" %lf ",a[j+rows*i]);

        }

        printf(" |\n");
    }
}

/* Multiplikation von b mit Q* */
void
qr_transform(int ld, double* qr, int m, int n, double* b) {



}

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in QR-Darstellung */
void
solve_qr_decomp(int ld, double* qr, int m, int n, double* b, double *x) {



}
