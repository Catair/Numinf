/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 		     	Serie 1 - Test mergesort		 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Oliver Heﬂ, Jan Witzany		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>	/* Fuer die Zufallseintraege  */

void
mergesort(int n, int* a){

  int m1, m2, i, k, l;
  int *b, *c;

  if(n > 1){
    m1 = (int) n/2;
    m2 = n - m1;

    /* Aufteilen der Liste a in b und c*/
    b = (int *) malloc(sizeof(int) * m1);
    c = (int *) malloc(sizeof(int) * m2);

    for(i=0; i<m1; i++){
        b[i] = a[i];
    }

    for(i=0; i<m2; i++){
        c[i] = a[i+m1];
    }

    mergesort(m1, b);
    mergesort(m2, c);

    i = 0; k = 0; l = 0;

    while(i < n){
        while(k < m1 && (l == m2 || b[k] <= c[l])){
            a[i] = b[k];
            k++;
            i++;
        }
        while(l < m2 && (k == m1 || c[l] <= b[k])){
            a[i] = c[l];
            l++;
            i++;
        }
    }

    /* Speicher wieder freigeben */
    free(c);
    free(b);
  }

}


/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */


int
main(void){

  int* a;	/* Vektor */
  int n;	/* Dimension des Problems */
  time_t t;	/* Zeitkonstante */
  int i;	/* Laufindex */

  time(&t);
  srand(t);	/* Zufallsgenerator */

  n = 10;

  /* ----------------------------------------
     Anlegen und befuellen Zufallsvektor
     ---------------------------------------- */

  a = (int*) malloc(sizeof(int) * n);
  for(i = 0; i < n; i++){
    a[i] = rand() % n;
  }

  printf(" Eingabevektor: \n");
  printf("( %d,", a[0]);
  for(i = 1; i < n-1; i++){
    printf(" %d,", a[i]);
  }
  printf(" %d)\n", a[n-1]);


  /* ----------------------------------------
     Anwenden mergesort
     ---------------------------------------- */

  printf("Starte mergesort\n");
  mergesort(n, a);

  printf(" Ausgabevektor: \n");
  printf("( %d,", a[0]);
  for(i = 1; i < n-1; i++){
    printf(" %d,", a[i]);
  }
  printf(" %d)\n", a[n-1]);

  /* Speicherfreigabe */
  free(a);

  return EXIT_SUCCESS;
}
