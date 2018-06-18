/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 	Serie 9 - Interpolation (Neville-Aitken-Verfahren /	 */
/*			  Tschebyscheff-interpolationspunkte)	 */
/* ------------------------------------------------------------- */
/*	Autoren:    	Oliver Heﬂ  & Jan Witzany	 		 */
/*	Versionsnummer: 					 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif

/* Der Parameter 'x' entspricht dem Punkt, an dem die Funktionen
   ausgewertet werden sollen, in der Variable 'data' koennen zusaetzliche
   Informationen zum Beispiel ein Parameter uebergeben werden.
   */
typedef double (*function)(double x, void *data);

double
sinus(double x, void *data);

double
f_test(double x, void *data);

double
f_test2(double x, void *data);

double
expo(double x, void *data);

/* P15
   Interpolationspolynom in y mit dem Neville-Aitken-Verfahren auswerten */
double
neville(double y, double* x, function f, void *fdata, int m);

/* Berechne aequidistante Interpolationspunkte */
void
aequidistant(double *x, int m, double a, double b);

/* P16
   Tschebyscheff-Interpolationspunkte im Intervall [a,b] berechnen */
void
tschebyscheff(double *x, int m, double a, double b);


/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int
main(void){

    int m, i, n;
    double y, a, b, h;
    double *xa, *xt;
    FILE *F;
    function f;
    double *param;
    double lambda;

    /* ----------------------------------------
    Startwerte und Speicher
    ---------------------------------------- */
    a = -5.0;
    b = 5.0;
    m = 10;
    y = 2.432;		/* Testpunkt */
    f = f_test;
    lambda = 50.0;
    param = &lambda;

    xa = (double*) calloc(m+1, sizeof(double));
    xt = (double*) calloc(m+1, sizeof(double));

    aequidistant(xa, m, a, b);
    tschebyscheff(xt, m, a, b);

    /* ------------------------------------------------------------
    *  Neville-Aitken-Verfahren
    * ------------------------------------------------------------ */
    printf("\n-----------Neville - aequidistant---------\n");
    printf("Interpolant in %.2f : %f\n", y, neville(y, xa, f, param, m));
    printf("Originalfunktion in %.2f : %f\n", y, f(y, param));
    printf("---------------------------\n");

    F = fopen("./solution_aequi.txt", "w");
    n = 100;
    h = (b-a)/n;
    for(i = 0; i < n+1; i++){
        fprintf(F, "%f %f\n", a + i*h, neville(a + i*h, xa, f, param, m));
    }
    fclose(F);

    printf("\n-----------Neville - Tschebycheff---------\n");
    printf("Interpolant in %.2f : %f\n", y, neville(y, xt, f, param, m));
    printf("Originalfunktion in %.2f : %f\n", y, f(y, param));
    printf("---------------------------\n");

    F = fopen("./solution_tscheb.txt", "w");
    n = 100;
    h = (b-a)/n;
    for(i = 0; i < n+1; i++){
        fprintf(F, "%f %f\n", a + i*h, neville(a + i*h, xt, f, param, m));
    }
    fclose(F);

    /* Speicherfreigabe */
    free(xa);
    free(xt);

    return EXIT_SUCCESS;
}

/* =============================== */
/*	  Hilfsfunktionen 	   */
/* =============================== */

double
neville(double y, double* x, function f, void *fdata, int m){
    int n,i,j;
    double erg;
    double *fj;

    fj = (double*) calloc(m+1, sizeof(double));

    for(j=0;j<=m;j++) {
        fj[j] = f(x[j], fdata);
    }

    for(n=1;n<=m;n++) {

        for(j=m;j>=n;j--) {
            i = j-n;
            fj[j] = ((y-x[i]) * fj[j] + (x[j]-y) * fj[j-1]) / (x[j]-x[i]);
        }
    }
    erg = fj[m];
    free(fj);

    return erg;
}

void
aequidistant(double *x, int m, double a, double b){
    int i;
    double dis, step;

    dis = b-a;
    step = dis/(m+1);

    x[0] = a;
    for(i=1;i<=m-1;i++) {
        x[i] = x[i-1] + step;
    }
    x[m] = b;
}

void
tschebyscheff(double *x, int m, double a, double b){
    int i;

    for(i=0;i<=m;i++) {
        x[i] = ((b+a)/2 + (b-a)/2) * cos(M_PI * ((2.0 * i + 1.0) / (2.0*m + 2.0)));
        printf("cos = %lf    ",cos(M_PI * ((2.0*i + 1.0) / (2.0*m + 2.0))));
        printf("tsche x %d = %lf\n",i,x[i]);
    }
}


/* =============================== */
/*   Testfunktionen 	   	   */
/* =============================== */


double
f_test(double x, void *data){
    double *lambda = data;
    return *lambda * exp(-1/2*x*x);
}

double f_test2(double x, void *data) {
    (void) data;
    return 1/(1+x*x);
}

double
sinus(double x, void *data){
    (void) data;
    return sin(x);
}

double
expo(double x, void *data) {
    (void) data;
    return exp(x);
}

