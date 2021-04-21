/*#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;

vfloat vec1, vec2 ;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 float res ;
 int i ;

 init_flop () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
  {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;
     
     start = _rdtsc () ;
        res = mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
  }
}*/

#include <stdio.h>
#include <x86intrin.h>

#ifndef MNBLAS_H
#define MNBLAS_H
#include "mnblas.h"
#endif

#ifndef COMPEXE_H
#define COMPEXE_H
#include "complexe.h"
#endif

#ifndef TESTUILS_H
#define TESTUILS_H
#include "testutils.h"
#endif

#define VECSIZE_RESULTAT    20
#define RAND_MAXIMUM        10

int main (int argc, char **argv) {
    float* V1s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float));
    double* V1d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t));

    float* V2s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float));
    double* V2d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double));
    complexe_float_t* V2c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* V2z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t));



    complexe_float_t tmpc;
    complexe_double_t tmpz;
    


    printf("                   TEST_DOT\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("---------- test 1 : \n");
    void_vector_sinit(V1s, 1,VECSIZE_RESULTAT);
    void_vector_sinit(V2s, 1,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_sdot : %f\n",mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1));
    printf("---------- test 2 : \n");
    void_vector_sinit(V2s, 2,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_sdot : %f\n",mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1));
    printf("---------- test 3 : \n");
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_sdot : %f\n",mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1));
    printf("---------- test 3 : \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_sdot : %f\n",mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1));

    printf("<------------------------------------------>\n                   double\n");

    printf("---------- test 1 : \n");
    void_vector_dinit(V1d, 1,VECSIZE_RESULTAT);
    void_vector_dinit(V2d, 1,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_ddot : %f\n",mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1));
    printf("---------- test 2 : \n");
    void_vector_dinit(V2d, 2,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_ddot : %f\n",mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1));
    printf("---------- test 3 : \n");
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_ddot : %f\n",mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1));
    printf("---------- test 4 : \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_ddot : %f\n",mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1));

    printf("<------------------------------------------>\n                   complexe_float_t (mncblas_cdotu_sub)\n");

    printf("---------- test 1 : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2 : \n");
    void_vector_cinit(V2c, gen_complexe_float(2,0),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2_bis : \n");
    void_vector_cinit(V2c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2_ter : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 3 : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 4 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotu_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);

    printf("<------------------------------------------>\n                   complexe_float_t (mncblas_cdotc_sub)\n");

    printf("---------- test 1 : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2 : \n");
    void_vector_cinit(V2c, gen_complexe_float(2,0),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2_bis : \n");
    void_vector_cinit(V2c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 2_ter : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(1,1),VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 3 : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);
    printf("---------- test 4 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    printf("mncblas_cdotc_sub : %.2f + %.2fi\n",tmpc.real,tmpc.imaginary);

    printf("<------------------------------------------>\n                   complexe_double_t (mncblas_zdotc_sub)\n");

    printf("---------- test 1 : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2 : \n");
    void_vector_zinit(V2z, gen_complexe_double(2,0),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2_bis : \n");
    void_vector_zinit(V2z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2_ter : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 3 : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 4 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotu_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);

    printf("<------------------------------------------>\n                   complexe_double_t (mncblas_zdotc_sub)\n");

    printf("---------- test 1 : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2 : \n");
    void_vector_zinit(V2z, gen_complexe_double(2,0),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2_bis : \n");
    void_vector_zinit(V2z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 2_ter : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(1,1),VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 3 : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);
    printf("---------- test 4 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    printf("mncblas_zdotc_sub : %.2f + %.2fi\n",tmpz.real,tmpz.imaginary);


    free(V1s);
    free(V1d);
    free(V1c);
    free(V1z);

    #define VECSIZE_FLOPS       100000
    #define NB_EXPE_VISIBLE     6
    #define NB_EXPE             1000
    #define NB_OPE_REEL         (2*VECSIZE_FLOPS) // PAR VECTEUR
    #define NB_OPE_COMPLEXE     (8*VECSIZE_FLOPS) // PAR VECTEUR
    #define NB_OPE_COMPLEXE_C   (9*VECSIZE_FLOPS) // PAR VECTEUR

    V1s = (float*)malloc(VECSIZE_FLOPS*sizeof(float));
    V1d = (double*)malloc(VECSIZE_FLOPS*sizeof(double));
    V1c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t));
    V2s = (float*)malloc(VECSIZE_FLOPS*sizeof(float));
    V2d = (double*)malloc(VECSIZE_FLOPS*sizeof(double));
    V2c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t));
    V2z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t));

    // V1s = vector_sinit(1,VECSIZE_FLOPS);
    // V1d = vector_dinit(1,VECSIZE_FLOPS);
    // V1c = vector_cinit(gen_complexe_float(1,1),VECSIZE_FLOPS);
    // V1z = vector_zinit(gen_complexe_double(1,1),VECSIZE_FLOPS);

    // V2s = vector_sinit(1,VECSIZE_FLOPS);
    // V2d = vector_dinit(1,VECSIZE_FLOPS);
    // V2c = vector_cinit(gen_complexe_float(1,1),VECSIZE_FLOPS);
    // V2z = vector_zinit(gen_complexe_double(1,1),VECSIZE_FLOPS);

    init_flop();

    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1);
        end = _rdtsc();
        calcul_flop("mncblas_sdot : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE (%d)\n",NB_EXPE);
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_sdot (VECSIZE_RESULTAT,V1s,1,V2s,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_sdot : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1);
        end = _rdtsc();
        calcul_flop("mncblas_ddot : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_ddot (VECSIZE_RESULTAT,V1d,1,V2d,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_ddot : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t (mncblas_cdotu_sub)\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
        end = _rdtsc();
        calcul_flop("mncblas_cdotu_sub : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t (mncblas_cdotu_sub) sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cdotu_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cdotu_sub : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t (mncblas_cdotc_sub)\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
        end = _rdtsc();
        calcul_flop("mncblas_cdotc_sub : ", NB_OPE_COMPLEXE_C ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t (mncblas_cdotc_sub) sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cdotc_sub(VECSIZE_RESULTAT,V1c,1,V2c,1,&tmpc);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cdotc_sub : ", NB_EXPE*NB_OPE_COMPLEXE_C ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t (mncblas_zdotu_sub)\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
        end = _rdtsc();
        calcul_flop("mncblas_zdotu_sub : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t (mncblas_zdotu_sub) sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zdotu_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zdotu_sub : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t (mncblas_zdotc_sub)\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
        end = _rdtsc();
        calcul_flop("mncblas_zdotc_sub : ", NB_OPE_COMPLEXE_C ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t (mncblas_zdotc_sub) sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zdotc_sub(VECSIZE_RESULTAT,V1z,1,V2z,1,&tmpz);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zdotc_sub : ", NB_EXPE*NB_OPE_COMPLEXE_C ,end-start);



    free(V1s);
    free(V1d);
    free(V1c);
    free(V1z);
}
