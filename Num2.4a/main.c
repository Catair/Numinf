/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 		     	Serie 2 - Test forward  		 */
/* ------------------------------------------------------------- */
/*	Autoren: 		    Oliver Heﬂ  Jan Witzany	 		 */
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

/* Loest ein lineares Gleichungssystem durch Vorwaertseinsetzen */
void
forward_subst(int ld, double* l, double* b, double *x){

    int i,j;

    for(j=0;j<ld;j++) {
        x[j] = b[j] / l[(ld+1)*j];
        for(i=(j+1);i<ld;i++) {
            b[i]= b[i]-l[i+j*ld]*x[j];
        }
    }
}

/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int main(void){

  int i, j;
  int n;
  double *a, *b, *x, *y;
  time_t t;
  time(&t);
  srand(t);

  n = 3;			/* Problemdimension*/

  /* ----------------------------------------
     Speicher abfordern
     ---------------------------------------- */
  a = (double *) malloc(sizeof(double) * n*n); /* Matrix */
  b = (double *) malloc(sizeof(double) * n);   /* Rechte Seite */
  x = (double *) malloc(sizeof(double) * n);   /* Loesungsvektor */
  y = (double *) malloc(sizeof(double) * n);   /* Pruefvektor */

  /* ----------------------------------------
     Befuellen mit Zufallseintraegen
     ---------------------------------------- */
  for(i = 0; i < n; i++){
    for(j = 0; j < i + 1; j++){
      set_entry(a, n, i, j, (rand() % 90)/10.0 + 1.0);  /* a_ii ungleich 0 */
    }
  }

  for(i=0;i<n;i++){
    set_entry(b,n,i,0,(rand() % 200)/20.0);
  }

  /* ----------------------------------------
     Ausgeben der Startwerte (Matrix / Vektor)
     ---------------------------------------- */
  printf("\n a:\n");
  print_matrix(a, n, n);

  printf("\n b:\n");
  print_matrix(b, n, 1);

  /* ----------------------------------------
     Loesen und ausgeben
     ---------------------------------------- */
  forward_subst(n, a, b, x);

  printf("\n x:\n");
  print_matrix(x, n, 1);

  /* ----------------------------------------
     Probe
     ---------------------------------------- */
  mvm(a, n, n, x, y);

  printf("\n y:\n");
  print_matrix(y, n, 1);

  /* Speicherfreigabe */
  free(a);
  free(b);
  free(x);
  free(y);

  return EXIT_SUCCESS;
}
