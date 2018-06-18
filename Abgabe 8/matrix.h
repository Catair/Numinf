/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 	Header File Matrix	 	 			 */
/* ------------------------------------------------------------- */
/*	Autoren: 		 				 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#ifndef MATRIX_H
#define MATRIX_H

/*----------------------------------------------------
 Matrix handeling
 -----------------------------------------------------*/

/* Eintrag (row,col) aus der Matrix a auslesen */
double
get_entry(double* a, int ld, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void
set_entry(double* a, int ld, int row, int col, double value);

/* Matrix Ausgabe (auf dem Bildschirm) */
void
print_matrix(double* a, int rows, int cols);

/*----------------------------------------------------
Creating matrices and arithmetic
------------------------------------------------------*/

/* Erstellen der Matrix des 1d-Modellproblems */
void
build_matrix1d(double* a, int n);

/* Matrix-Vektor-Multiplikation (Ax = y) */
void
mvm(double* a, int rows, int cols, double* x, double* y);

/* Matrix-Vektor-Multiplikation Tridiagonalmatrix */
void
mvm_tridiag(int trans, double alpha, int n, double *l, double *d, double *u,
	    double *x, double *y);

double
skaldxminusy(double* x, double* y, double delta, int n);

/* Berechnet die euklidische Norm */
double
norm_2(double* x, int n);

/* Berechnet das euklidische Skalarprodukt */
double
scalar_2(double* x, double* y, int n);

/* Vektoriteration auf a anwenden und bei Genauigkeit eps abbrechen */
void
power_adaptive(double* a, int n, double* x, double eps);

/* Inverse Iteration auf a anwenden und bei Genauigkeit eps abbrechen */
void
invit_adaptive(double* a, int n, double* x, double eps);

/* Die Arrays 'l,d,u' entsprechenden den Neben- und der Haupdiagonalen der Marix.
   Der Parameter 'shift' enthaelt den Shift der Matrix.
   Die Groesse 'steps' wird genutzt, um die Anzahl der Iterationsschritte fertzulegen.
   Im Pointer 'eigenvalue' wird der errechnete Eigenwert gespeichert.
   Der Pointer 'res' wird genutzt, um die Norm des Residuums || Ax - lambda x|| nach aussen zu tragen. */

void
inverse_iteration_withshift_tridiag(int n, double *l, double *d, double *u, double *x, double shift, int steps, double *eigenvalue, double *res);

/* Die Arrays 'l,d,u' entsprechenden den Neben- und der Haupdiagonalen der Marix.
   Der Parameter 'shift' enthaelt den Shift der Matrix.
   Die Groesse 'steps' wird genutzt, um die Anzahl der Iterationsschritte fertzulegen.
   Im Pointer 'eigenvalue' wird der errechnete Eigenwert gespeichert.
   Der Pointer 'res' wird genutzt, um die Norm des Residuums || Ax - lambda x|| nach aussen zu tragen. */

void
rayleigh_iteration_tridiag(int n, double *l, double *d, double *u, double *x, double shift, int steps, double *eigenvalue, double *res);

/*----------------------------------------------------
LR-decomposition and solving triangular systems
------------------------------------------------------*/

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void
backward_subst(int ld, double* r, double* b, double* x);

void
forward_substlr(int ld, double* l, double* b, double *x);

/* Loesen mit pivotiserter LR-Zerlegung */
void
solve_lr_pivot(int ld, const double* a, int* p, double* b, double* x);

/* LR-Zerlegung mit Pivot */
void
lr_pivot(int ld, double* a, int* p);

/* LR-Zerlegung */
void
lr_decomp(int ld, double* a);

/* Loesen mit LR-Zerlegung */
void
solve_lr_decomp(int ld, const double* a, const double* b, double* x);

/* LR-Zerlegung Tridiagonalmatrix */
void
lr_decomp_tridiag(int n, double *l, double *d, double *u);

/* Loesen mit LR-Zerlegung Tridiagonalmatrix */
void
lr_solve_tridiag(int n, double *l, double *d, double *u, double *b, double *x);

/*----------------------------------------------------
QR-decomposition and solving
------------------------------------------------------*/

/* Berechnet die QR-Zerlegung der (rows x cols)-Matrix a */
void
qr_decomp(int ld, double* a, int rows, int cols);

/* Multiplikation von b mit Q* */
void
qr_transform(int ld, double* qr, int m, int n, double* b);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in QR-Darstellung */
void
solve_qr_decomp(int ld, double* qr, int m, int n, double* b, double *x);

#endif
