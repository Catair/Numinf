/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 	Serie 6 - Test Newton- und Heron-Verfahren 		 */
/* ------------------------------------------------------------- */
/*	Autoren:  				 		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------
 Funktionsprototypen f :R -> R, df `: R->R
------------------------------------------------------------*/

/* Der Parameter 'x' entspricht dem Punkt, an dem die Funktionen
   asugewertet werden soll, in der Variable 'data' koennen zusaetzliche
   Informationen zum Beispiel ein Parameter uebergeben werden.
   */

typedef double (*function)(double x, void *data);

/* Funktionen aus Aufgabe T11 (a)*/
/* Beispielfunktion f in x auswerten */
double
f1(double x, void *data);

/* Ableitung von f in x auswerten */
double
df1(double x, void *data);

/* Weitere Beispielfunktion (frei waehlbar) */
double
f2(double x, void *data);

double
df2(double x, void *data);

/* Newton-Verfahren mit Genauigkeit eps */
double
newton(function f, void *fdata, function df, void *dfdata, double x0, double eps);

/* Heron-Verfahren mit Genauigkeit eps */
double
heron(double x0, double eps);

/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */

int
main(void){

    //Newton Test mit T11
    printf("f(x) = x^3 - 3x^2 + x -3\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert 0.0 ist:\n");
    printf("%lf\n\n",newton(f1,NULL,df1,NULL,0.0,1e10));

    printf("\nf(x) = x^3 - 3x^2 + x -3\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert 0.5 ist:\n");
    printf("%lf\n\n",newton(f1,NULL,df1,NULL,0.5,1e10));

    printf("\nf(x) = x^3 - 3x^2 + x -3\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert 4.0 ist:\n");
    printf("%lf\n\n",newton(f1,NULL,df1,NULL,4.0,1e10));

    //Newton Test mit eigener Fkt
    printf("\nf(x) = sin(x)\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert 0.0 ist:\n");
    printf("%lf\n\n",newton(f2,NULL,df2,NULL,0.0,1e10));

    printf("\nf(x) = sin(x)\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert -45.0 ist:\n");
    printf("%lf\n\n",newton(f2,NULL,df2,NULL,-45.0,1e10));

    printf("\nf(x) = sin(x)\n");
    printf("Die Näherung der Nullstelle von f(x) mit Startwert 110.0 ist:\n");
    printf("%lf\n\n",newton(f2,NULL,df2,NULL,110.0,1e10));

    //Heron Aufruf
    printf("\nDie Wurzel von w ist: ");
    printf("%lf",heron(1.0,0.0001));

    return EXIT_SUCCESS;
}

/* =============================== */
/*	  Hilfsfunktionen 	   */
/* =============================== */

double
f1(double x, void *data){

    return (x*x*x)-3*(x*x)+x-3;

}

double
df1(double x, void *data){

    return 3*(x*x)-6*x+1;

}

double
f2(double x, void *data){

    return sin(x);

}

double
df2(double x, void *data){

    return cos(x);

}

double
newton(function f, void *fdata, function df, void *dfdata, double x0, double eps){

    double gr = 1+eps;
    double x;

    while(gr >= eps) {
        gr = f(x, fdata)/df(x, dfdata);
        x = x - gr;
    }

    return x;
}

double
heron(double x0, double eps){

    int m = 0;
    double x = x0;
    const double w = 25.0;

    while(fabs(x*x - w) >= eps) {

        x = (x +(w/x))/2.0;
        m++;
    }

    return x;
}
