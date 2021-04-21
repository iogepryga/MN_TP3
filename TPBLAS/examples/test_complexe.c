/*#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{
  complexe_float_t c1= {1.0, 2.0} ;
  complexe_float_t c2= {3.0, 6.0} ;
  complexe_float_t c3 ;

  complexe_double_t cd1 ;
  complexe_double_t cd2 ;
  complexe_double_t cd3 ;

  unsigned long long int start, end ;
  int i ;

  init_flop () ;
  
  c3 = add_complexe_float (c1, c2) ;

  printf ("c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

  cd1 = (complexe_double_t) {10.0, 7.0} ;
  cd2 = (complexe_double_t) {25.0, 32.0} ;

  cd3 = add_complexe_double (cd1, cd2) ;

  printf ("cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

  start =_rdtsc () ;
  
  for (i = 0 ; i < NB_FOIS; i++)
    {
      cd3 = add_complexe_double (cd1, cd2) ;
    }

  end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;
  calcul_flop ("calcul complexe ", NB_FOIS*2, end-start) ;


  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe float\n----------------------------------------------\n");
  c3 = mult_complexe_float(c1,c2);
  printf("c1 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c1);
  printf("c2 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c1,c1);
  printf("c1 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c2);
  printf("c2 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++) {
      c3 = mult_complexe_float (c1, c2) ;
  }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe float", NB_FOIS*6, end-start) ;
  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe double\n----------------------------------------------\n");
  cd3 = mult_complexe_double(cd1,cd2);
  printf("cd1 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd1);
  printf("cd2 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd1,cd1);
  printf("cd1 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd2);
  printf("cd2 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++){
    cd3 = mult_complexe_double (cd1, cd2) ;
  }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe double", NB_FOIS*6, end-start) ;
  exit (0) ;
}*/

#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#define NB_FOIS       4194304
#define RAND_MAXIMUM  10

#include "flop.h"


complexe_float_t gen_complexe_float (const float real, const float imaginary) {
    complexe_float_t tmp; tmp.real = real; tmp.imaginary = imaginary;
    return tmp;
}
complexe_double_t gen_complexe_double (const double real, const double imaginary) {
    complexe_double_t tmp; tmp.real = real; tmp.imaginary = imaginary;
    return tmp;
}

float rand_float(const int max) {
    return ((float)rand()/(float)(RAND_MAX)) * max;
}

double rand_double(const int max) {
    return ((float)rand()/(float)(RAND_MAX)) * max;
}

complexe_float_t rand_complexe_float_t(const int max) {
    complexe_float_t tmp;
    tmp.real = rand_float(max);
    tmp.imaginary = rand_float(max);
    return tmp;
}

complexe_double_t rand_complexe_double_t(const int max) {
    complexe_double_t tmp;
    tmp.real = rand_double(max);
    tmp.imaginary = rand_double(max);
    return tmp;
}

int main (int argc, char **argv)
{
  srand(_rdtsc());
 complexe_float_t c1 = gen_complexe_float(1,2);
 complexe_float_t c2 = gen_complexe_float(3,6);
 complexe_float_t c3;

 complexe_double_t cd1 = gen_complexe_double(10,7);
 complexe_double_t cd2 = gen_complexe_double(25,32);
 complexe_double_t cd3;

 unsigned long long int start, end, sum;
 register int i ;

 init_flop () ;
 

 printf ("c1 = %.2f + %.2fi\n", c1.real, c1.imaginary);
 printf ("c2 = %.2f + %.2fi\n", c2.real, c2.imaginary);
 c3 = add_complexe_float (c1, c2) ;
 printf ("c3 = c1 + c2 = %.2f + %.2fi\n", c3.real, c3.imaginary);


 printf ("cd1 = %.2f + %.2fi\n", cd1.real, cd1.imaginary);
 printf ("cd2 = %.2f + %.2fi\n", cd2.real, cd2.imaginary);
 cd3 = add_complexe_double (cd1, cd2) ;
 printf ("cd3 = cd1 + cd2 = %.2f + %.2fi\n", cd3.real, cd3.imaginary);

 /*start =_rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd3 = add_complexe_double (cd1, cd2) ;
   }

 end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;
  calcul_flop ("calcul complexe ", NB_FOIS*2, end-start) ;*/

  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe float\n----------------------------------------------\n");
  c3 = mult_complexe_float(c1,c2);
  printf("c1 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c1);
  printf("c2 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c1,c1);
  printf("c1 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c2);
  printf("c2 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  sum = 0;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      c1 = rand_complexe_float_t(RAND_MAXIMUM);
      c2 = rand_complexe_float_t(RAND_MAXIMUM);
      start =_rdtsc () ;
      mult_complexe_float (c1, c2) ;
      end = _rdtsc () ;
      sum += end-start;
    }
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe float", NB_FOIS*6, sum) ;
  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe double\n----------------------------------------------\n");
  cd3 = mult_complexe_double(cd1,cd2);
  printf("cd1 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd1);
  printf("cd2 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd1,cd1);
  printf("cd1 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd2);
  printf("cd2 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  sum = 0;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      cd1 = rand_complexe_double_t(RAND_MAXIMUM);
      cd2 = rand_complexe_double_t(RAND_MAXIMUM);
      start =_rdtsc () ;
      mult_complexe_double (cd1, cd2) ;
      end = _rdtsc () ;
      sum += end-start;
    }
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe double", NB_FOIS*6, sum) ;

  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Division complexe float\n----------------------------------------------\n");
  c3 = div_complexe_float(c1,c2);
  printf("c1 / c2 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c2,c1);
  printf("c2 / c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c1,c1);
  printf("c1 / c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c2,c2);
  printf("c2 / c2 = %f+%fi\n",c3.real,c3.imaginary);
  sum = 0;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      c1 = rand_complexe_float_t(RAND_MAXIMUM);
      c2 = rand_complexe_float_t(RAND_MAXIMUM);
      start =_rdtsc () ;
      div_complexe_float (c1, c2) ;
      end = _rdtsc () ;
      sum += end-start;
    }
  printf("Calcul flop : \n");
  calcul_flop ("div complexe float", NB_FOIS*6, sum) ;
  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Division complexe double\n----------------------------------------------\n");
  cd3 = div_complexe_double(cd1,cd2);
  printf("cd1 / cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd2,cd1);
  printf("cd2 / cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd1,cd1);
  printf("cd1 / cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd2,cd2);
  printf("cd2 / cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  sum = 0;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      cd1 = rand_complexe_double_t(RAND_MAXIMUM);
      cd2 = rand_complexe_double_t(RAND_MAXIMUM);
      start =_rdtsc () ;
      div_complexe_double (cd1, cd2) ;
      end = _rdtsc () ;
      sum += end-start;
    }
  printf("Calcul flop : \n");
  calcul_flop ("div complexe double", NB_FOIS*6, sum) ;
  exit (0) ;
}
