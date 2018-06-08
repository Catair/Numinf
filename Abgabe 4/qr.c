/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/*	Serie 4 - Test QR-Zerlegung (Givens)			 */
/* ------------------------------------------------------------- */
/*	Autoren: 				 		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define min(a,b)  (((a) < (b)) ? (a) : (b))

/* Eintrag (row,col) aus der Matrix a auslesen */
double
get_entry(double* a, int ld, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void
set_entry(double* a, int ld, int row, int col, double value);

/* Matrix Ausgabe (auf dem Bildschirm) */
void
print_matrix(double* a, int rows, int cols);

/* Berechnet die QR-Zerlegung der (rows x cols)-Matrix a */
void
qr_decomp(int ld, double* a, int rows, int cols);

/* Testen der QR-Zerlegung */
void
test_qr(int ld, double* a, double* qr, int rows, int cols);


/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int main(void){

  int i, j, rows, cols, min;
  double *a, *qr;
  time_t t;
  time(&t);
  srand(t);

  rows = 3;		/* Matrixzeilen */
  cols = 4;		/* Matrixspalten */

  /* ----------------------------------------
     Speicher anfordern und Befuellen
     ---------------------------------------- */
  a = (double *) malloc(sizeof(double)* rows * cols);
  qr = (double *) malloc(sizeof(double) * rows * cols );

  for(i = 0; i < rows; i++){
    for(j = 0; j < cols; j++){
      set_entry(a, rows, i, j, (rand() % 200)/20.0);
      set_entry(qr, rows, i, j, get_entry(a, rows, i, j));
    }
  }

  /* Ausgabe Informationen und Matrix */
  printf("\n rows = %d, cols = %d \n",rows, cols);
  min = (rows > cols ? cols : rows);
  printf("Minimum = %d \n", min);
  printf("\n a:\n");
  print_matrix(qr, rows, cols);

  /* ----------------------------------------
     QR - Zerlegung durchfuehren und testen
     ---------------------------------------- */

  qr_decomp(min(rows,cols), qr, rows, cols);

  printf("\n qr:\n");
  print_matrix(qr, rows, cols);

  test_qr(min(rows,cols), a, qr, rows, cols);

  /* Speicherfreigabe */
  free(a);
  free(qr);

  return EXIT_SUCCESS;
}


/* =============================== */
/*	  Hilfsfunktionen 	   */
/* =============================== */


double
get_entry(double* a, int ld, int row, int col){

    return a[row+ld*col];

}

void
set_entry(double* a, int ld, int row, int col, double value){

    a[row+ld*col] = value;

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

void
qr_decomp(int ld, double* a, int rows, int cols){

    int i, j, k;
    double c, s, alpha, delta, tau;

    for(k=0;k < min(rows,cols);k++) {

        for(i=k+1;i<rows;i++) {

            if(get_entry(a, ld, i,k) == 0.0) {

                delta = 1;
                c = 1;
                s = 0;

            }else{

                if(fabs(get_entry(a, ld, k, k)) >= fabs(get_entry(a, ld, i, k))) {

                    tau = get_entry(a, ld, i, k) / get_entry(a, ld, k, k);
                    delta = tau / sqrt(pow(tau,2) + 1);
                    s = delta;
                    c = sqrt(1 - pow(s,2));

                }else{

                tau = get_entry(a, ld, k, k) / get_entry(a, ld, i, k);
                delta = sqrt(pow(tau,2))/tau;
                c = 1 / delta;
                s = sqrt(1 - pow(c,2));

                }
            }

            set_entry(a, ld, k, k, ((c * get_entry(a, ld, k, k)) + (s * get_entry(a, ld, i, k))));

            for(j=k+1;j<cols;j++) {

                alpha = get_entry(a, ld, k, j);
                set_entry(a, ld, k, j, (c*alpha + s*get_entry(a, ld, i, j)));
                set_entry(a, ld, i, j, ((-1)*s*alpha + c*get_entry(a, ld, i, j)));
            }
        }
    }
}

void
test_qr(int ld, double* a, double* qr, int rows, int cols){

    int i,j;
    static int RUND = 100000;
    double b, r;

    for(i=0;i<cols;i++) {

        a = 0;
        r = 0;

        for(j=0;j<rows;j++) {          //berechnen der eukl. Norm in zwei Schritten 1. Summe unter der Wurzel

            b = b + pow(get_entry(a, ld, j, i), 2);
            r = r + pow(get_entry(qr, ld, j, i), 2);
        }

        //Wurzel der Normberechnung
        b = sqrt(b);
        r = sqrt(r);

        //Vergleich der Normen auf n Nachkommastellen
        b = b * RUND;
        r = r * RUND;

        if((int)b != (int)r) {
            printf("Fehler in Spalte %d Normen nicht identisch",i);
        }
    }
}
