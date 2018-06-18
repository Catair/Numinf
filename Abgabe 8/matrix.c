/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 	C File Matrix	 	 				 */
/* ------------------------------------------------------------- */
/*	Autoren: 		 				 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"

#define min(a,b)  (((a) < (b)) ? (a) : (b))

/*----------------------------------------------------
 Matrix handeling
 -----------------------------------------------------*/

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

/*----------------------------------------------------
Creating matrices and arithmetic
------------------------------------------------------*/

void
build_matrix1d(double* a, int n){

    int i,j;

    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {

            if(i==j) {

                set_entry(a, n, i, j, 2.0);

            }else{

                if((i-j == 1) || (j-i == 1)) {

                    set_entry(a, n, i, j, -1.0);
                }else{

                    set_entry(a, n, i, j, 0.0);
                }
            }
        }
    }
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

void
mvm_tridiag(int trans, double alpha, int n, double *l, double *d, double *u,
	    double *x, double *y){

  int i;

  if(n == 0){
    return;
  }

  if(n == 1) {
    y[0] += alpha * d[0] * x[0];
    return;
  }

  if(trans) {
    y[0] += alpha * (d[0] * x[0] + l[0] * x[1]);
    for(i = 1; i < n-1; i++){
      y[i] += alpha * (u[i-1] * x[i-1] + d[i] * x[i] + l[i] * x[i+1]);
    }
    y[i] += alpha * (u[i-1] * x[i-1] + d[i] * x[i]);
  }
  else {
    y[0] += alpha * (d[0] * x[0] + u[0] * x[1]);
    for(i = 1; i < n-1; i++){
      y[i] += alpha * (l[i-1] * x[i-1] + d[i] * x[i] + u[i] * x[i+1]);
    }
    y[i] += alpha * (l[i-1] * x[i-1] + d[i] * x[i]);
  }
}

double
skaldxminusy(double* x, double* y, double delta, int n) {

    int i;
    double *erg;
    erg = (double *) malloc(sizeof(double) * n);

    for(i=0;i<n;i++) {

        erg[i] = (delta*x[i]) - y[i];
    }

    return norm_2(erg, n);
}

double
norm_2(double* x, int n){

    int i;
    double s = 0.0;

    for(i=0;i<n;i++) {
        s = s + (x[i]*x[i]);
    }
    return sqrt(s);
}


double
scalar_2(double* x, double* y, int n){

    int i;
    double s = 0.0;

    for(i=0;i<n;i++) {

        s = s + (x[i]*y[i]);
    }

    return s;
}

void
power_adaptive(double* a, int n, double* x, double eps){

    if(norm_2(x,n) != 0.0) {    //Die Norm ist genau dann 0.0 wenn alle Elemente von x auch 0.0 sind

        int i;
        double gamma, delta, *y;

        y = (double *) malloc(sizeof(double) * n);

        gamma = norm_2(x, n);

        for(i=0;i<n;i++) {
            x[i] = x[i]/gamma;
        }

        mvm(a, n, n, x, y);

        delta = scalar_2(x, y, n);

        while(skaldxminusy(x, y, delta, n) > (eps*norm_2(y, n))) {

            gamma = norm_2(y, n);

            for(i=0;i<n;i++){
                   x[i] = y[i]/gamma;
            }

            mvm(a, n, n, x, y);
            delta =  scalar_2(x, y, n);
        }

        free(y);

    }else{
        printf("X ist Nullvektor!");
    }
}


void
invit_adaptive(double* a, int n, double* x, double eps){

    if(norm_2(x,n) != 0.0) {    //Die Norm ist genau dann 0.0 wenn alle Elemente von x auch 0.0 sind
        int i,j;
        double gamma, lambda, *y, *xs, *b;

        y = (double *) malloc(sizeof(double) * n);
        xs = (double *) malloc(sizeof(double) * n);
        b = (double *) malloc(sizeof(double) * n*n);

        //Sicherungskopien von a und x für qr Zerlegung
        for(i=0;i<n;i++){
            for(j=0;j<n;j++) {

                set_entry(b, n, i, j, get_entry(a, n, i, j));
            }
            xs[i] = x[i];
        }

        gamma = norm_2(x, n);

        for(i=0;i<n;i++) {
            x[i] = x[i]/gamma;
        }

        qr_decomp(n, b, n, n);
        solve_qr_decomp(n, b, n, n, xs, y);

        lambda = scalar_2(x,y,n);

        while(skaldxminusy(x, y, lambda, n) > (eps*norm_2(y, n))) {

            gamma = norm_2(y, n);

            for(i=0;i<n;i++){
                x[i] = y[i]/gamma;
                xs[i]=x[i];
            }

            solve_qr_decomp(n, b, n, n, xs, y);

            lambda =  scalar_2(x, y, n);
        }

        free(b);
        free(xs);
        free(y);

    }else{
        printf("X ist Nullvektor!");
    }
}


void
inverse_iteration_withshift_tridiag(int n, double *l, double *d, double *u, double *x, double shift, int steps, double *eigenvalue, double *res){

    int i,j;
    double gamma, *y;

    y = (double *) malloc(sizeof(double) * n);

    gamma = norm_2(x, n);

    for(i=0;i<n;i++) {
        x[i] = x[i]/gamma;
    }

    //Shift
    for(i=0;i<n;i++) {
        d[i] = d[i] - shift;
    }

    lr_decomp_tridiag(n, l, d, u);
    lr_solve_tridiag(n, l, d, u, x, y);

    *eigenvalue = scalar_2(x,y,n);

    for(j=0;j<steps;j++) {

        gamma = norm_2(y, n);

        for(i=0;i<n;i++){
            x[i] = y[i]/gamma;
         }

        lr_solve_tridiag(n, l, d, u, x, y);

        *eigenvalue =  scalar_2(x, y, n);
    }

    //y wird hier wiederverwendet und überschrieben
    mvm_tridiag(0, 1, n, l, d, u, x, y);

    for(i=0;i<n;i++) {
        y[i] = y[i] - x[i]* *eigenvalue;
    }

    *res = norm_2(y,n);

    free(y);
}

void
rayleigh_iteration_tridiag(int n, double *l, double *d, double *u, double *x, double shift, int steps, double *eigenvalue, double *res){

    int i,j,t=0;
    double gamma, *y;

    y = (double *) malloc(sizeof(double) * n);

    gamma = norm_2(x, n);

    for(i=0;i<n;i++) {
        x[i] = x[i]/gamma;
    }

    //Shift
    for(i=0;i<n;i++) {
        d[i] = d[i] - shift;
    }

    lr_decomp_tridiag(n, l, d, u);
    lr_solve_tridiag(n, l, d, u, x, y);

    *eigenvalue = scalar_2(x,y,n);

    shift = (1 / *eigenvalue) + shift;

    for(j=0;j<steps;j++) {

        gamma = norm_2(y, n);

        for(i=0;i<n;i++){
            x[i] = y[i]/gamma;
        }

        lr_solve_tridiag(n, l, d, u, x, y);

        *eigenvalue =  scalar_2(x, y, n);
        shift = (1 / *eigenvalue) + shift;
    }

  //y wird hier wiederverwendet und überschrieben
    mvm_tridiag(0, 1, n, l, d, u, x, y);

    for(i=0;i<n;i++) {
        y[i] = y[i] - x[i]* *eigenvalue;
    }

    *res = norm_2(y,n);

    printf("test %d\n",t);
    t++;

    free(y);

    printf("test %d\n",t);
    t++;
}

/*----------------------------------------------------
LR-decomposition and solving triangular systems
------------------------------------------------------*/


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
forward_substlr(int ld, double* l, double* b, double *x){

    int i,j;

    for(j=0;j<ld;j++) {
        x[j] = b[j] / l[(ld+1)*j];
        for(i=(j+1);i<ld;i++) {
            b[i]= b[i]-l[i+j*ld]*x[j];
        }
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

void
lr_pivot(int ld, double* a, int* p){
    int i,is,j,k;       //is ist eine Hilfsvariable für die Zeilenvertauschung
    double g;

    for(k=0;k<ld;k++){

        is = k;

        for(i=(k+1);i<ld;i++) {         //Finde maximales Element

            if(fabs(get_entry(a, ld, i, k)) > get_entry(a, ld, is, k)) {
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

void
solve_lr_pivot(int ld, const double* a, int* p, double* b, double* x){
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
                    set_entry(l, ld, i, j, get_entry((double*)a, ld, i, j));
                }else{
                    set_entry(l, ld, i, j, 0.0);
                }
            }
        }
    }

    /* Ly = b wird ausgewertet*/
    forward_substlr(ld, l, b, y);

    backward_subst(ld, (double*)a, y, x );

    free(l);
    free(y);
}

/* Spezialfall Tridiagonal */

void
lr_decomp_tridiag(int n, double *l, double *d, double *u){

    int i;

    for(i=0;i<n-1;i++) {

        l[i] = l[i]/d[i];
        d[i+1] = d[i+1] - l[i]*u[i];
    }
}

void
lr_solve_tridiag(int n, double *l, double *d, double *u, double *b, double *x){

    int i;

    //modifizierter Forward sub Algorithmus
    for(i=0;i<n-1;i++) {
        x[i] = b[i];
        b[i+1] = b[i+1]-l[i]*x[i];
    }
    x[n]=b[n];

    //Vorbereitung für backward Sub
    for(i=0;i<n;i++) {
        b[i] = x[i];
    }

    //modifizierter Backward sub Alg.
    for(i=n;i>0;i--) {
        x[i] = b[i]/d[i];
        b[i-1] = b[i-1] - u[i-1]*x[i];
    }
}


/*----------------------------------------------------
QR-decomposition and solving
------------------------------------------------------*/

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


void
qr_transform(int ld, double* qr, int m, int n, double* b){

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
solve_qr_decomp(int ld, double* qr, int m, int n, double* b, double *x){

    qr_transform(ld, qr, m, n, b);
    backward_subst(ld, qr, b, x);
}
