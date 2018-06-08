/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 		     	Serie 2 - Test Matrix			 */
/* ------------------------------------------------------------- */
/*	Autoren: 			Oliver Heﬂ  Jan Witzany			 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Eintrag (row,col) aus der Matrix a auslesen */

double
get_entry(double* a, int ld, int row, int col){

    return a[row+ld*col];

}

/* Eintrag (row,col) der Matrix a als value setzen */
/* Vektor a hat die Zellen 0 bis 8 */
void
set_entry(double* a, int ld, int row, int col, double value){

    a[row+ld*col] = value;

}

/* Matrix Ausgabe (auf dem Bildschirm) */

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


/* Matrix-Vektor-Multiplikation (Ax = y) */
// rows = n     cols = m
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



/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int
main(void){

  int MODULO = 3;
  int n, m;			/* Matrixgroessen */
  int i, j;			/* Lauf Indices */
  double *a, *x, *y;		/* Matrizen und Vektoren */
  time_t t;
  time(&t);
  srand(t);

  n = 3;			/* Zeilen */
  m = 6;			/* Spalten */

  /* ----------------------------------------
     Speicher abfordern
     ---------------------------------------- */
  a = (double *) malloc(sizeof(double) *(n * m));
  x = (double *) malloc(sizeof(double) * m);
  y = (double *) malloc(sizeof(double) * n);


  /* ----------------------------------------
     Befuellen mit Zufallseintraegen
     ---------------------------------------- */
  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      set_entry(a, n, i, j, (rand() % 200)/20.0);
    }
  }

  for(i = 0; i < m; i++){
    set_entry(x, m, i, 0, (rand() % 200)/20.0);
  }

  /* ----------------------------------------
     Ausgeben der Startwerte (Matrix / Vektor)
     ---------------------------------------- */

  printf("\n a:\n");
  print_matrix(a, n, m);

  printf("\n x:\n");
  print_matrix(x, m, 1);

  /* ----------------------------------------
     Matrix-Vektor-Multiplikation
     ---------------------------------------- */

  mvm(a, n, m, x, y);

  printf("\ny:\n");
  print_matrix(y, n, 1);

  /* Speicherfreigabe */
  printf("Speicherfreigabe\n");
  free(a);
  free(x);
  free(y);

  return EXIT_SUCCESS;
}

