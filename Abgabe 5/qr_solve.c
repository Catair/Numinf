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
#define	M_PI		3.14159265358979323846

void
test_qr(double* a, double* x, double* bt, double* b, int rows, int cols);

void
test_matrix1d(double* a, double* x, double* b,int square);

//Baut u''(x) für Test
void
testhilf_2(double* x, int n);

//Baut u(x) für Test
void
testhilf(double* x, int n);

/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int main(void){

/* =============================== */
/*	  Aufgabe 9b 	   */
/* =============================== */

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

/* ---------------------------------------
QR - Zerlegung durchfuehren und testen
---------------------------------------- */

    qr_decomp(min(rows,cols), qr, rows, cols);

    printf("\n QR Zerlegung:\n");
    print_matrix(qr, rows, cols);

    solve_qr_decomp(rows, qr, rows, cols, b, x);

    printf("\n x:\n");
    print_matrix(x, rows, 1);

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

    test_qr(qr, x, bt, b, rows, cols);

    printf("\n Testfall qr*x = bt\n");
    printf(" Matrix a:\n");
    print_matrix(qr, rows, cols);
    printf("\n Testvektor bt:\n");
    print_matrix(bt, rows, 1);
    printf("\n    QR Test ende\n\n\n");

/* Speicherfreigabe */
    free(b);
    free(bt);
    free(x);
    free(qr);

/* =============================== */
/*	  Aufgabe 9c 	   */
/* =============================== */
    int square;
    double *a, *y, *c;

    square = 8;     /* Quadratmatrix */

    a = (double *) malloc(sizeof(double) * square * square);
    y = (double *) malloc(sizeof(double) * square);
    c = (double *) malloc(sizeof(double) * square);

    build_matrix1d(a,square);

    print_matrix(a,square,square);

    printf("\nStart Test.\n");

    test_matrix1d(a, y, c, square);

    free(a);
    free(y);
    free(c);

    return EXIT_SUCCESS;
}

/*----------------------------------------------------
Tests
------------------------------------------------------*/

void
test_qr(double* a, double* x, double* bt, double *b, int rows, int cols) {

    int i;
    double eps = 1e-8;

    mvm(a, rows, cols, x, bt);

    for(i=0;i<rows;i++) {

        if(fabs(b[i])-fabs(bt[i]) > eps) {
            printf("\nFehler in Zeile %d\n",i);
            printf("|%lf| - |%lf|",b[i],bt[i]);
        }
    }
}

void
test_matrix1d(double* a, double* x, double* b,int square) {

    int i;

    double eps=1e-8;
    double *xabs;   //analytische Lsg für x;

    xabs = (double *) malloc(sizeof(double) * square);

    printf("Testfunktion u(x)= exp(-pi^2t)*sin(pi*x)");
    printf("\n\nTestmarixgroesse: %d",square);

    testhilf_2(b,square);

    qr_decomp(square, a, square, square);
    solve_qr_decomp(square, a, square, square, b, x);
    testhilf(xabs,square);

    for(i=0;i<square;i++) {
        if(fabs(xabs[i])-fabs(x[i]) > eps) {

            printf("\nFehler in Zeile %d\n",i);
        }
    }
}

//Baut u(x) für Test
void
testhilf(double* x, int n) {

    int i;
    double h = 1.0 / (n+1.0), xnow = h;

    for(i=0;i<n;i++) {

        x[i]= (exp(-1.0*M_PI*M_PI) * sin(M_PI*xnow));
        xnow = xnow + h;
    }
}

//Baut u''(x) für Test
void
testhilf_2(double* x, int n) {

    int i=0;
    double h = 1.0 / (n + 1.0), xnow = h;

        x[0]= M_PI*M_PI * exp(-1.0*M_PI*M_PI) * -1.0 * sin(M_PI*xnow);

        for(i=1;i<n;i++) {

        x[i]= M_PI*M_PI * exp(-1.0*M_PI*M_PI) * -1.0 * sin(M_PI*xnow);
        xnow = xnow + h;
    }
}


