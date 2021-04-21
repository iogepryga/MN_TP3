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

#define M_RESULTAT 3
#define N_RESULTAT 4
#define K_RESULTAT 2
#define RAND_MAXIMUM 10

int main (int argc, char **argv) {
    float* M1s = (float*)malloc(M_RESULTAT*N_RESULTAT*sizeof(float));
    float* M2s = (float*)malloc(N_RESULTAT*K_RESULTAT*sizeof(float));
    float* M3s = (float*)malloc(M_RESULTAT*K_RESULTAT*sizeof(float));
    float* M3s_copy = (float*)malloc(M_RESULTAT*K_RESULTAT*sizeof(float));
    double* M1d = (double*)malloc(M_RESULTAT*N_RESULTAT*sizeof(double));
    double* M2d = (double*)malloc(N_RESULTAT*K_RESULTAT*sizeof(double));
    double* M3d = (double*)malloc(M_RESULTAT*K_RESULTAT*sizeof(double));
    double* M3d_copy = (double*)malloc(M_RESULTAT*K_RESULTAT*sizeof(double));
    complexe_float_t* M1c = (complexe_float_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_float_t));
    complexe_float_t* M2c = (complexe_float_t*)malloc(N_RESULTAT*K_RESULTAT*sizeof(complexe_float_t));
    complexe_float_t* M3c = (complexe_float_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_float_t));
    complexe_float_t* M3c_copy = (complexe_float_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* M1z = (complexe_double_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_double_t));
    complexe_double_t* M2z = (complexe_double_t*)malloc(N_RESULTAT*K_RESULTAT*sizeof(complexe_double_t));
    complexe_double_t* M3z = (complexe_double_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_double_t));
    complexe_double_t* M3z_copy = (complexe_double_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_double_t));
    
    
    
    complexe_float_t tmpc = gen_complexe_float(2,0);
    complexe_double_t tmpz = gen_complexe_double(2,0);

    


    printf("                   TEST_GEMM\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("-------------------------- test 1 (M3s=2*M1s*M2s+2*M3s) (MNCblasRowMajor): \n");
    void_matrix_sinit(M1s, 1,M_RESULTAT,N_RESULTAT);
    void_matrix_sinit(M2s, 1,N_RESULTAT,K_RESULTAT);
    void_matrix_sinit(M3s, 1,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1s : "); matrix_print(M1s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2s : "); matrix_print(M2s,TYPE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3s : "); matrix_print(M3s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_sgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,2,M1s,1,M2s,1,2,M3s,1);
    printf("---- Apres :\n");
    printf("M3s : "); matrix_print(M3s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("-------------------------- test 3 (M3s=2*M1s*M2s+2*M3s) (MNCblasRowMajor): \n");
    void_matrix_sinit_rand(M1s, RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_matrix_sinit_rand(M2s, RAND_MAXIMUM,N_RESULTAT,K_RESULTAT);
    void_matrix_sinit_rand(M3s, RAND_MAXIMUM,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1s : "); matrix_print(M1s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2s : "); matrix_print(M2s,TYPE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3s : "); matrix_print(M3s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_sgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,2,M1s,1,M2s,1,2,M3s,1);
    printf("---- Apres :\n");
    printf("M3s : "); matrix_print(M3s,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("<------------------------------------------>\n                   double\n");

    printf("-------------------------- test 1 (M3d=2*M1d*M2d+2*M3d) (MNCblasRowMajor): \n");
    void_matrix_dinit(M1d, 1,M_RESULTAT,N_RESULTAT);
    void_matrix_dinit(M2d, 1,N_RESULTAT,K_RESULTAT);
    void_matrix_dinit(M3d, 1,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1d : "); matrix_print(M1d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2d : "); matrix_print(M2d,TYPE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3d : "); matrix_print(M3d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_dgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,2,M1d,1,M2d,1,2,M3d,1);
    printf("---- Apres :\n");
    printf("M3d : "); matrix_print(M3d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("-------------------------- test 3 (M3d=2*M1d*M2d+2*M3d) (MNCblasRowMajor): \n");
    void_matrix_dinit_rand(M1d, RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_matrix_dinit_rand(M2d, RAND_MAXIMUM,N_RESULTAT,K_RESULTAT);
    void_matrix_dinit_rand(M3d, RAND_MAXIMUM,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1d : "); matrix_print(M1d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2d : "); matrix_print(M2d,TYPE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3d : "); matrix_print(M3d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_dgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,2,M1d,1,M2d,1,2,M3d,1);
    printf("---- Apres :\n");
    printf("M3d : "); matrix_print(M3d,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("-------------------------- test 1 (M3c=2*M1c*M2c+2*M3c) (MNCblasRowMajor): \n");
    void_matrix_cinit(M1c, gen_complexe_float(1,0),M_RESULTAT,N_RESULTAT);
    void_matrix_cinit(M2c, gen_complexe_float(1,0),N_RESULTAT,K_RESULTAT);
    void_matrix_cinit(M3c, gen_complexe_float(1,0),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1c : "); matrix_print(M1c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2c : "); matrix_print(M2c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpc,M1c,1,M2c,1,&tmpc,M3c,1);
    printf("---- Apres :\n");
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    printf("-------------------------- test 1_bis (M3c=2*M1c*M2c+2*M3c) (MNCblasRowMajor): \n");
    void_matrix_cinit(M1c, gen_complexe_float(0,1),M_RESULTAT,N_RESULTAT);
    void_matrix_cinit(M2c, gen_complexe_float(0,1),N_RESULTAT,K_RESULTAT);
    void_matrix_cinit(M3c, gen_complexe_float(0,1),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1c : "); matrix_print(M1c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2c : "); matrix_print(M2c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpc,M1c,1,M2c,1,&tmpc,M3c,1);
    printf("---- Apres :\n");
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    printf("-------------------------- test 1_ter (M3c=2*M1c*M2c+2*M3c) (MNCblasRowMajor): \n");
    void_matrix_cinit(M1c, gen_complexe_float(1,1),M_RESULTAT,N_RESULTAT);
    void_matrix_cinit(M2c, gen_complexe_float(1,1),N_RESULTAT,K_RESULTAT);
    void_matrix_cinit(M3c, gen_complexe_float(1,1),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1c : "); matrix_print(M1c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2c : "); matrix_print(M2c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpc,M1c,1,M2c,1,&tmpc,M3c,1);
    printf("---- Apres :\n");
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("-------------------------- test 3 (M3c=2*M1c*M2c+2*M3c) (MNCblasRowMajor): \n");
    void_matrix_cinit_rand(M1c, RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_matrix_cinit_rand(M2c, RAND_MAXIMUM,N_RESULTAT,K_RESULTAT);
    void_matrix_cinit_rand(M3c, RAND_MAXIMUM,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1c : "); matrix_print(M1c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2c : "); matrix_print(M2c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpc,M1c,1,M2c,1,&tmpc,M3c,1);
    printf("---- Apres :\n");
    printf("M3c : "); matrix_print(M3c,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("-------------------------- test 1 (M3z=2*M1z*M2z+2*M3z) (MNCblasRowMajor): \n");
    void_matrix_zinit(M1z, gen_complexe_double(1,0),M_RESULTAT,N_RESULTAT);
    void_matrix_zinit(M2z, gen_complexe_double(1,0),N_RESULTAT,K_RESULTAT);
    void_matrix_zinit(M3z, gen_complexe_double(1,0),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1z : "); matrix_print(M1z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2z : "); matrix_print(M2z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpz,M1z,1,M2z,1,&tmpz,M3z,1);
    printf("---- Apres :\n");
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    printf("-------------------------- test 1_bis (M3z=2*M1z*M2z+2*M3z) (MNCblasRowMajor): \n");
    void_matrix_zinit(M1z, gen_complexe_double(0,1),M_RESULTAT,N_RESULTAT);
    void_matrix_zinit(M2z, gen_complexe_double(0,1),N_RESULTAT,K_RESULTAT);
    void_matrix_zinit(M3z, gen_complexe_double(0,1),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1z : "); matrix_print(M1z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2z : "); matrix_print(M2z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpz,M1z,1,M2z,1,&tmpz,M3z,1);
    printf("---- Apres :\n");
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    printf("-------------------------- test 1_ter (M3z=2*M1z*M2z+2*M3z) (MNCblasRowMajor): \n");
    void_matrix_zinit(M1z, gen_complexe_double(1,1),M_RESULTAT,N_RESULTAT);
    void_matrix_zinit(M2z, gen_complexe_double(1,1),N_RESULTAT,K_RESULTAT);
    void_matrix_zinit(M3z, gen_complexe_double(1,1),M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1z : "); matrix_print(M1z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2z : "); matrix_print(M2z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpz,M1z,1,M2z,1,&tmpz,M3z,1);
    printf("---- Apres :\n");
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);

    printf("-------------------------- test 3 (M3z=2*M1z*M2z+2*M3z) (MNCblasRowMajor): \n");
    void_matrix_zinit_rand(M1z, RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_matrix_zinit_rand(M2z, RAND_MAXIMUM,N_RESULTAT,K_RESULTAT);
    void_matrix_zinit_rand(M3z, RAND_MAXIMUM,M_RESULTAT,K_RESULTAT);
    printf("---- Avant :\n");
    printf("M1z : "); matrix_print(M1z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("M2z : "); matrix_print(M2z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,N_RESULTAT,K_RESULTAT);
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);
    mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,K_RESULTAT,&tmpz,M1z,1,M2z,1,&tmpz,M3z,1);
    printf("---- Apres :\n");
    printf("M3z : "); matrix_print(M3z,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,K_RESULTAT);






    free(M1s);
    free(M2s);
    free(M3s);
    free(M3s_copy);
    free(M1d);
    free(M2d);
    free(M3d);
    free(M3d_copy);
    free(M1c);
    free(M2c);
    free(M3c);
    free(M3c_copy);
    free(M1z);
    free(M2z);
    free(M3z);
    free(M3z_copy);

    #define M_FLOPS         100
    #define N_FLOPS         10
    #define K_FLOPS         10
    #define NB_EXPE         1000
    #define NB_OPE_REEL     (M_FLOPS*N_FLOPS*K_FLOPS*2+M_FLOPS*3)
    #define NB_OPE_COMPLEXE (M_FLOPS*N_FLOPS*K_FLOPS*8+M_FLOPS*14)


    M1s = (float*)malloc(M_RESULTAT*N_RESULTAT*sizeof(float));
    M2s = (float*)malloc(N_RESULTAT*K_RESULTAT*sizeof(float));
    M3s = (float*)malloc(M_RESULTAT*K_RESULTAT*sizeof(float));
    M1d = (double*)malloc(M_RESULTAT*N_RESULTAT*sizeof(double));
    M2d = (double*)malloc(N_RESULTAT*K_RESULTAT*sizeof(double));
    M3d = (double*)malloc(M_RESULTAT*K_RESULTAT*sizeof(double));
    M1c = (complexe_float_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_float_t));
    M2c = (complexe_float_t*)malloc(N_RESULTAT*K_RESULTAT*sizeof(complexe_float_t));
    M3c = (complexe_float_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_float_t));
    M1z = (complexe_double_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_double_t));
    M2z = (complexe_double_t*)malloc(N_RESULTAT*K_RESULTAT*sizeof(complexe_double_t));
    M3z = (complexe_double_t*)malloc(M_RESULTAT*K_RESULTAT*sizeof(complexe_double_t));

    init_flop();

    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n");
    unsigned long long int start, end ; 
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_sgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_FLOPS,N_FLOPS,K_FLOPS,2,M1s,1,M2s,1,2,M3s,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_sgemm : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_dgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_FLOPS,N_FLOPS,K_FLOPS,2,M1d,1,M2d,1,2,M3d,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_dgemm : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_FLOPS,N_FLOPS,K_FLOPS,&tmpc,M1c,1,M2c,1,&tmpc,M3c,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cgemm : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M_FLOPS,N_FLOPS,K_FLOPS,&tmpz,M1z,1,M2z,1,&tmpz,M3z,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zgemm : ", NB_EXPE*NB_OPE_COMPLEXE ,end-start);


    // free(M1s);
    // free(M2s);
    // free(M3s);
    // free(M1d);
    // free(M2d);
    // free(M3d);
    // free(M1c);
    // free(M2c);
    // free(M3c);
    // free(M1z);
    // free(M2z);
    // free(M3z);
}