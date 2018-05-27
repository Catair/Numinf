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

#define min(a,b)  (((a) < (b)) ? (a) : (b))

// Test
void
test_qr(int ld, double* qr, double* x, double* bt, int rows, int cols);

/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int main(void){
  int rows, cols, min;
  double *x, *b, *bt, *qr;

  rows = 3;		/* Matrixzeilen */
  cols = 3;		/* Matrixspalten */

  /* ----------------------------------------
     Speicher anfordern und Befuellen
     ---------------------------------------- */
  qr = (double *) malloc(sizeof(double) * rows * cols );
  b = (double *) malloc(sizeof(double) * rows);
  bt = (double *) malloc(sizeof(double) * rows);
  x = (double *) malloc(sizeof(double) * rows);

  set_entry(qr, rows, 0, 0, 3.0);
  set_entry(qr, rows, 0, 1, -9.0);
  set_entry(qr, rows, 0, 2, 7.0);
  set_entry(qr, rows, 1, 0, -4.0);
  set_entry(qr, rows, 1, 1, -13.0);
  set_entry(qr, rows, 1, 2, -1.0);
  set_entry(qr, rows, 2, 0, 0.0);
  set_entry(qr, rows, 2, 1, -20.0);
  set_entry(qr, rows, 2, 2, -35.0);

  b[0] = -3.0;
  b[1] = -21.0;
  b[2] = -20.0;

  /* Ausgabe Informationen und Matrix */
  printf("\nrows = %d, cols = %d \n",rows, cols);
  min = (rows > cols ? cols : rows);
  printf("Minimum = %d \n", min);
  printf("\n a:\n");
  print_matrix(qr, rows, cols);
  printf("\n b:\n");
  print_matrix(b, rows, 1);

  /* ----------------------------------------
     QR - Zerlegung durchfuehren und testen
     ---------------------------------------- */

  qr_decomp(min(rows,cols), qr, rows, cols);

  printf("\n QR Zerlegung:\n");
  print_matrix(qr, rows, cols);

  solve_qr_decomp(rows, qr, rows, cols, b, x);

  printf("\n x:\n");
  print_matrix(x, rows, 1);

  mvm(qr, rows, cols, x, bt);

  printf("\n Testfall qr*x = bt\n bt:\n");
  print_matrix(bt, rows, 1);

  /* Speicherfreigabe */
  free(b);
  free(bt);
  free(x);
  free(qr);


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

void
mvm(double* a, int rows, int cols, double* x, double* y){

    int i,j;

       for(i=0;i<rows;i++) {       //y wird 0 gesetzt

        y[i]=0;

    }

    for(i=0;i<rows;i++) {

        for(j=0;j<cols;j++) {
            y[i]=y[i]+a[i+j*rows]*x[j];
        }
    }
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

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in QR-Darstellung */
void
solve_qr_decomp(int ld, double* qr, int m, int n, double* b, double *x) {

    qr_transform(ld, qr, m, n, b);
    backward_subst(ld, qr, b, x);
}

/* Multiplikation von b mit Q* */
void
qr_transform(int ld, double* qr, int m, int n, double* b) {

    int i,k;
    double delta, alpha, c, s;

    for(k = 0;k < min(m,n);k++) {

        for(i = (k+1);i < m;i++) {

            delta = get_entry(qr, ld, i, k);

            if(delta == 1.0) {

                c = 1.0;
                s = 0.0;

            }else{
                if(fabs(delta) < 1.0) {

                    s = delta;
                    c = sqrt(1.0 - (s*s));

                }else{

                    c = 1.0/delta;
                    s = sqrt(1.0 - (c*c));

                }
            }

            alpha = b[k];
            b[k] = c*alpha + s*b[i];
            b[i] = (-1.0)*s*alpha + c*b[i];
        }
    }
}

void
qr_decomp(int ld, double* a, int rows, int cols){

    int i, j, k;
    double c, s, alpha, delta, tau;

    for(k=0;k < min(rows,cols);k++) {

        for(i=k+1;i<rows;i++) {

            if(get_entry(a, ld, i,k) == 0.0) {

                delta = 1.0;
                c = 1.0;
                s = 0.0;

            }else{

                if(fabs(get_entry(a, ld, k, k)) >= fabs(get_entry(a, ld, i, k))) {

                    tau = get_entry(a, ld, i, k) / get_entry(a, ld, k, k);
                    delta = tau / sqrt(tau*tau + 1.0);
                    s = delta;
                    c = sqrt(1.0 - s*s);

                }else{

                tau = get_entry(a, ld, k, k) / get_entry(a, ld, i, k);
                delta = sqrt((tau*tau)+1.0)/tau;
                c = 1.0 / delta;
                s = sqrt(1.0 - c*c);

                }
            }

            set_entry(a, ld, k, k, ((c * get_entry(a, ld, k, k)) + (s * get_entry(a, ld, i, k))));
            set_entry(a, ld, i, k, delta);

            for(j=k+1;j<cols;j++) {

                alpha = get_entry(a, ld, k, j);
                set_entry(a, ld, k, j, (c*alpha + s*get_entry(a, ld, i, j)));
                set_entry(a, ld, i, j, ((-1)*s*alpha + c*get_entry(a, ld, i, j)));
            }
        }
    }
}

// Test
void
test_qr(int ld, double* qr, double* x, double* bt, int rows, int cols) {

    int i, j;
    double temp = 0.0;

    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
                temp = temp + (get_entry(qr, ld, i, j) * x[j]);
        }
        bt[i] = temp;
    }
}
