/*---------------------------------------------------------------*/
/*	Numerische Mathematik in den Ingenieurwissenschaften 	 */
/* 		     	Serie 1 - Test gnomesort		 */
/* ------------------------------------------------------------- */
/*	Autoren: 						 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>	/* Fuer die Zufallseintraege  */

#define max(a,b)  (((a) > (b)) ? (a) : (b))

int gnomesort(int n, int* a){

  int op=0;
  int i=1;
  int h;

  printf("\n");

  while(i < n) {

     if(a[i-1] <= a[i]){
        i=i+1;

        op++;
    }else{

        h = a[i-1];
        a[i-1] = a[i];
        a[i] = h;

        i = max(1, i-1);

        op++;
        }
    }
    printf("\n");
    return op;

}


/* =============================== */
/*	  Hauptprogramm 	   */
/* =============================== */


int main(void){

  int* a;	/* Vektor */
  int n;	/* Dimension des Problems */
  time_t t;	/* Zeitkonstante */
  int i;	/* Laufindex */
  int op;	/* Operationen */

  time(&t);
  srand(t);	/* Zufallsgenerator */

  n = 4;

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
     Anwenden gnomesort
     ---------------------------------------- */

  printf("Starte gnomesort\n");
  op = gnomesort(n, a);
  printf("Es wurden %d Operationen benoetigt\n", op);

  printf(" Ausgabevektor: \n");
  printf("( %d,", a[0]);
  for(i = 1; i < n-1; i++){
    printf(" %d,", a[i]);
  }
  printf(" %d)\n", a[n-1]);


  /* ----------------------------------------
     Worst-Case
     ---------------------------------------- */

  printf("Worst-Case\n");
  for(i = 0; i < n; i++){
    a[i] = n - i;
  }
  printf(" Eingabevektor: \n");
  printf("( %d,", a[0]);
  for(i = 1; i < n-1; i++){
    printf(" %d,", a[i]);
  }
  printf(" %d)\n", a[n-1]);

  printf("Starte gnomesort\n");
  op = gnomesort(n, a);
  printf("Es wurden %d Operationen benoetigt\n", op);

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
