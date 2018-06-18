/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 		     	Serie 3 - Test LR-Zerlegung 		 */
/* ------------------------------------------------------------- */
/*	Autoren: 		Oliver Heﬂ Jan Witzany	 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double
get_entry(double* a, int ld, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void
set_entry(double* a, int ld, int row, int col, double value);

/* Matrix Ausgabe (auf dem Bildschirm) */
void
print_matrix(double* a, int rows, int cols);

/* Berechnet die LR-Zerlegung der Matrix a und ueberspeichert a */
void
lr_decomp(int ld, double* a);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in LR-Darstellung ohne Pivotsuche */
void
solve_lr_decomp(int ld, const double* a, const double* b, double* x);

/* Loest ein lineares Gleichungssystem durch Vorwaertseinsetzen */
void
forward_substlr(int ld, double* l, double* b, double *x);

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void
backward_subst(int ld, double* r, double* b, double* x);

/* Testet die Korrektheit des Vektors x */
void
test_gls(int ld, double* a, double* b, double* x);

/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int
main(void){

  int i;
  int n;
  double *a, *b, *x;

  n = 3;					/* Problemdimension*/

  /* ----------------------------------------
     Speicher anfordern und befuellen
     ---------------------------------------- */
  a = (double *) malloc(sizeof(double) * n*n); /* Matrix */
  b = (double *) malloc(sizeof(double) * n);   /* Rechte Seite */
  x = (double *) malloc(sizeof(double) * n);   /* Loesungsvektor */

  set_entry(a, n, 0, 0, 3.0);
  set_entry(a, n, 0, 1, 2.0);
  set_entry(a, n, 0, 2, 1.0);
  set_entry(a, n, 1, 0, 6.0);
  set_entry(a, n, 1, 1, 6.0);
  set_entry(a, n, 1, 2, 3.0);
  set_entry(a, n, 2, 0, 9.0);
  set_entry(a, n, 2, 1, 10.0);
  set_entry(a, n, 2, 2, 6.0);

  b[0] = 3.0;
  b[1] = 6.0;
  b[2] = 11.0;

  /* Matrix ausgeben */
  printf("\n a:\n");
  print_matrix(a,n,n);


  /* ----------------------------------------
     LR-Zerlegung und Loesen
     ---------------------------------------- */
  lr_decomp(n, a);

  /* Matrix ausgeben */
  printf("\n LR-Darstellung:\n");
  print_matrix(a, n, n);

  /* Rechte Seite ausgeben */
  printf("\n b:\n");
  for(i = 0; i < n; i++){
	printf("%.2f\n", b[i]);
  }

  /* Gleichungssystem loesen */
  solve_lr_decomp(n, a, b, x);

  /* Loesung ausgeben */
  printf("\n x:\n");
  for(i = 0; i < n; i++){
	printf("%.2f\n", x[i]);
  }

  /*Wiederherstellung von A und b f¸r den Test, ob x Ax=b lˆst*/
  set_entry(a, n, 0, 0, 3.0);
  set_entry(a, n, 0, 1, 2.0);
  set_entry(a, n, 0, 2, 1.0);
  set_entry(a, n, 1, 0, 6.0);
  set_entry(a, n, 1, 1, 6.0);
  set_entry(a, n, 1, 2, 3.0);
  set_entry(a, n, 2, 0, 9.0);
  set_entry(a, n, 2, 1, 10.0);
  set_entry(a, n, 2, 2, 6.0);

  b[0] = 3.0;
  b[1] = 6.0;
  b[2] = 11.0;

  test_gls(n, a, b, x);

  /* Speicherfreigabe */
  free(a);
  free(b);
  free(x);

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
lr_decomp(int ld, double* a){

    int i,j,k;

    for(k=0;k<ld;k++){

        for(i=(k+1);i<ld;i++) {
            set_entry(a, ld, i, k, (get_entry(a, ld, i, k)/get_entry(a, ld, k, k)));
        }

        for(i=(k+1);i<ld;i++) {

            for(j=(k+1);j<ld;j++) {
                set_entry(a, ld, i, j, (get_entry(a, ld, i, j)-(get_entry(a, ld, i, k) * get_entry(a, ld, k, j))));
            }
        }
    }
}

void
solve_lr_decomp(int ld, const double* a, const double* b, double* x){
    int i,j;
    double *l, *y, *bs;

    l = (double *) calloc(ld*ld, sizeof(double));
    y = (double *) calloc(ld, sizeof(double));
    bs = (double *) calloc(ld, sizeof(double));

    /* Teilmatrix L wird erzeugt */
    for(i=0;i<ld;i++) {
        for(j=0;j<ld;j++) {

            if(i==j) {
                set_entry(l, ld, i, j, 1.0);
            }else{
                if(i>j) {
                    set_entry(l, ld, i, j, get_entry((double*)a, ld, i, j));
                }else{
                    set_entry(l, ld, i, j, 0.0);
                }
            }
        }
        bs[i] = b[i];
    }

    /* Ly = b wird ausgewertet*/
    forward_substlr(ld, l, bs, y);

    backward_subst(ld, (double*)a, y, x );

    free(bs);
    free(l);
    free(y);
}

/* Loest ein lineares Gleichungssystem durch Vorwaertseinsetzen */
void
forward_substlr(int ld, double* l, double* b, double *x){

    int i,j;

    for(j=0;j<ld;j++) {
        x[j] = b[j] / l[(ld+1)*j];
        for(i=(j+1);i<ld;i++) {
            b[i]= b[i]-l[i+j*ld]*x[j];
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

/* Testet die Korrektheit des Vektors x */
void
test_gls(int ld, double* a, double* b, double* x) {
    int i,j;
    double* bs;  //b solved Ergebnis des Gleichungssystems Ax

    bs = (double *) calloc(ld, sizeof(double));

    for(i=0;i<ld;i++) {

        for(j=0;j<ld;j++) {

            bs[i] = bs[i] + (get_entry(a, ld, i, j)*x[j]);
        }
    }

    for(i=0;i<ld;i++) {

        if(bs[i] != b[i]) {
            printf("Fehler im Algorithmus!\n");
            print_matrix(b, ld, 1);
            printf("\n");
            print_matrix(bs, ld, 1);
            break;
        }
    }

    free(bs);
}
