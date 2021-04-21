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
    float* V1s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float)); float* V2s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float));
    double* V1d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double)); double* V2d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t)); complexe_float_t* V2c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t)); complexe_double_t* V2z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t));
    complexe_float_t tmpc; tmpc.real = 2; tmpc.imaginary = 0;
    complexe_double_t tmpz; tmpz.real = 2; tmpz.imaginary = 0;

    


    printf("                   TEST_AXPY\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("---------- test 1 (V2s=2*V1s+V2s): \n");
    void_vector_sinit(V1s, 1,VECSIZE_RESULTAT);
    void_vector_sinit(V2s, 2,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    mnblas_saxpy(VECSIZE_RESULTAT,2,V1s,1,V2s,1);
    printf("---- Apres :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2s=2*V1s+V2s): \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    V1s[VECSIZE_RESULTAT-1] = 0;
    V2s[VECSIZE_RESULTAT-1] = 1;
    printf("---- Avant :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    mnblas_saxpy(VECSIZE_RESULTAT,2,V1s,1,V2s,1);
    printf("---- Apres :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   double\n");

    printf("---------- test 1 (V2s=2*V1s+V2s): \n");
    void_vector_dinit(V1d, 1,VECSIZE_RESULTAT);
    void_vector_dinit(V2d, 2,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_daxpy(VECSIZE_RESULTAT,2,V1d,1,V2d,1);
    printf("---- Apres :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2d=2*V1d+V2d): \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    V1s[VECSIZE_RESULTAT-1] = 0;
    V2s[VECSIZE_RESULTAT-1] = 1;
    printf("---- Avant :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_daxpy(VECSIZE_RESULTAT,2,V1d,1,V2d,1);
    printf("---- Apres :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("---------- test 1 (V2d=2*V1d+V2d): \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit2(V2c,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mnblas_caxpy(VECSIZE_RESULTAT,&tmpc,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2d=2*V1d+V2d): \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1c[VECSIZE_RESULTAT-1] = gen_complexe_float(0,0);
    V2c[VECSIZE_RESULTAT-1] = gen_complexe_float(1,0);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mnblas_caxpy(VECSIZE_RESULTAT,&tmpc,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("---------- test 1 (V2d=2*V1d+V2d): \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit2(V2z,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_zaxpy(VECSIZE_RESULTAT,&tmpz,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2d=2*V1d+V2d): \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1z[VECSIZE_RESULTAT-1] = gen_complexe_double(0,0);
    V2z[VECSIZE_RESULTAT-1] = gen_complexe_double(1,0);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_zaxpy(VECSIZE_RESULTAT,&tmpz,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);


    free(V1s);
    free(V1d);
    free(V1c);
    free(V1z);
    free(V2s);
    free(V2d);
    free(V2c);
    free(V2z);

    #define VECSIZE_FLOPS       100000
    #define NB_EXPE_VISIBLE     6
    #define NB_EXPE             1000
    #define NB_OPE_REEL         (2*VECSIZE_FLOPS)
    #define NB_OPE_COMPLEXE     (8*VECSIZE_FLOPS)

    V1s = (float*)malloc(VECSIZE_FLOPS*sizeof(float)), V2s = (float*)malloc(VECSIZE_FLOPS*sizeof(float));
    V1d = (double*)malloc(VECSIZE_FLOPS*sizeof(double)), V2d = (double*)malloc(VECSIZE_FLOPS*sizeof(double));
    V1c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t)), V2c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t)), V2z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t));

    init_flop();

    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_saxpy(VECSIZE_FLOPS,2,V1s,1,V2s,1);
        end = _rdtsc();
        calcul_flop("mnblas_saxpy : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_saxpy(VECSIZE_FLOPS,2,V1s,1,V2s,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_saxpy : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_daxpy(VECSIZE_FLOPS,2,V1d,1,V2d,1);
        end = _rdtsc();
        calcul_flop("mnblas_daxpy : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_daxpy(VECSIZE_FLOPS,2,V1d,1,V2d,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_daxpy : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_caxpy(VECSIZE_FLOPS,&tmpc,V1c,1,V2c,1);
        end = _rdtsc();
        calcul_flop("mnblas_caxpy : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_caxpy(VECSIZE_FLOPS,&tmpc,V1c,1,V2c,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_caxpy : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_zaxpy(VECSIZE_FLOPS,&tmpz,V1z,1,V2z,1);
        end = _rdtsc();
        calcul_flop("mnblas_zaxpy : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_zaxpy(VECSIZE_FLOPS,&tmpz,V1z,1,V2z,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_zaxpy : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);



    free(V1s);
    free(V1d);
    free(V1c);
    free(V1z);
    free(V2s);
    free(V2d);
    free(V2c);
    free(V2z);
}