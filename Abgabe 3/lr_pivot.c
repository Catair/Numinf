/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/*	Serie 3 - Test LR-Zerlegung mit Pivotsuche 		 */
/* ------------------------------------------------------------- */
/*	Autoren: 		Oliver Heﬂ Jan Witzany	 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double
get_entry(double* a, int ld, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void
set_entry(double* a, int ld, int row, int col, double value);

/* Matrix Ausgabe (auf dem Bildschirm) */
void
print_matrix(double* a, int rows, int cols);

/* Berechnet die LR-Zerlegung der Matrix a mit Pivotsuche, ueberspeichert a und merkt sich die Vertauschungen in p*/
void
lr_pivot(int ld, double* a, int* p);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in LR-Darstellung mit Pivotsuche */
void
solve_lr_pivot(int ld, const double* a, int* p, double* b, double* x);

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

  int n;
  int i;
  double *a, *b, *x;
  int *p;

  n = 3;					/* Problemdimension*/

  /* ----------------------------------------
     Speicher anfordern und befuellen
     ---------------------------------------- */
  a = (double *) malloc(sizeof(double) * n*n);/* Matrix */
  b = (double *) malloc(sizeof(double) * n);  /* Rechte Seite */
  x = (double *) malloc(sizeof(double) * n);  /* Loesungsvektor */
  p = (int *) malloc(sizeof(int) * n);        /* Permutationsvektor */

  set_entry(a, n, 0, 0, 0.0);
  set_entry(a, n, 0, 1, 1.0);
  set_entry(a, n, 0, 2, 3.0);
  set_entry(a, n, 1, 0, 4.0);
  set_entry(a, n, 1, 1, 6.0);
  set_entry(a, n, 1, 2, 0.0);
  set_entry(a, n, 2, 0, 2.0);
  set_entry(a, n, 2, 1, 6.0);
  set_entry(a, n, 2, 2, 3.0);

  b[0] = 2.5;
  b[1] = -2.0;
  b[2] = 0.5;

  /* Matrix ausgeben */
  printf("\n a:\n");
  print_matrix(a, n, n);

  /* ----------------------------------------
     LR-Zerlegung und Loesen
     ---------------------------------------- */
  lr_pivot(n, a, p);

  /* Matrix ausgeben */
  printf("\n LR-Darstellung:\n");
  print_matrix(a, n, n);

  /* Permutationsvektor ausgeben */
  printf("\n p:\n");
  for(i = 0; i < n; i++){
	printf(" %d \n", p[i]);
  }

  /* Rechte Seite ausgeben */
  printf("\n b:\n");
  for(i = 0; i < n; i++){
	 printf("%.2f\n", b[i]);
  }

  /* Gleichungssystem loesen */
  solve_lr_pivot(n, a, p, b, x);

  /* Loesung ausgeben */
  printf("\n x:\n");
  for(i = 0; i < n; i++){
    printf("%.2f\n", x[i]);
  }

  /*Wiederherstellung von A und b f¸r den Test, ob x Ax=b lˆst*/
  set_entry(a, n, 0, 0, 0.0);
  set_entry(a, n, 0, 1, 1.0);
  set_entry(a, n, 0, 2, 3.0);
  set_entry(a, n, 1, 0, 4.0);
  set_entry(a, n, 1, 1, 6.0);
  set_entry(a, n, 1, 2, 0.0);
  set_entry(a, n, 2, 0, 2.0);
  set_entry(a, n, 2, 1, 6.0);
  set_entry(a, n, 2, 2, 3.0);

  b[0] = 2.5;
  b[1] = -2.0;
  b[2] = 0.5;

  test_gls(n, a, b, x);

  /* Speicherfreigabe */
  free(a);
  free(b);
  free(x);
  free(p);

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



void lr_pivot(int ld, double* a, int* p){

    int i,is,j,k;       //is ist eine Hilfsvariable f¸r die Zeilenvertauschung
    double g;

    for(k=0;k<ld;k++){

        is = k;

        for(i=(k+1);i<ld;i++) {         //Finde maximales Element

            if(fabs(get_entry(a, ld, i, k)) > (fabs(get_entry(a, ld, is, k)))) {
                is = i;
            }
        }
        p[k] = is;

        for(j=0;j<ld;j++) {             //Vertausche Zeilen

            g = get_entry(a, ld, k, j);
            set_entry(a, ld, k, j, get_entry(a, ld, is, j));
            set_entry(a, ld, is, j, g);

        }

        for(i=(k+1);i<ld;i++) {         // LR Zerlegung

            set_entry(a, ld, i, k, (get_entry(a, ld, i, k)/get_entry(a, ld, k, k)));
        }

        for(i=(k+1);i<ld;i++) {

            for(j=(k+1);j<ld;j++) {
                set_entry(a, ld, i, j, (get_entry(a, ld, i, j)-(get_entry(a, ld, i, k) * get_entry(a, ld, k, j))));
            }
        }
    }
}

void solve_lr_pivot(int ld, const double* a, int* p, double* b, double* x){
    int i,j;
    double *l, *y, g;

    l = (double *) calloc(ld*ld, sizeof(double));
    y = (double *) calloc(ld, sizeof(double));

    /* b wird an p angepasst */
    for(i=0;i<ld;i++) {
        if(i!=p[i]) {

            g = b[i];
            b[i] = b[p[i]];
            b[p[i]] = g;

        }
    }

    printf("\n pb:\n");
    for(i = 0; i < ld; i++){

        printf("%.2f\n", b[i]);
    }

    /* Teilmatrix L wird erzeugt */
    for(i=0;i<ld;i++) {
        for(j=0;j<ld;j++) {

            if(i==j) {
                set_entry(l, ld, i, j, 1.0);
            }else{
                if(i>j) {
                    set_entry(l, ld, i, j, get_entry(a, ld, i, j));
                }else{
                    set_entry(l, ld, i, j, 0.0);
                }
            }
        }
    }

    /* Ly = b wird ausgewertet*/
    forward_substlr(ld, l, b, y);

    for(i=0;i<ld;i++) {
        for(j=0;j<ld;j++) {

            if(i>j) {
                set_entry(a, ld, i, j, 0.0);
            }
        }
    }

    backward_subst(ld, a, y, x );

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
