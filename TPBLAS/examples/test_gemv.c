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

#define M_RESULTAT 2
#define N_RESULTAT 4
#define RAND_MAXIMUM 10

int main (int argc, char **argv) {
    int MAX_SIZE = (M_RESULTAT < N_RESULTAT ? N_RESULTAT : M_RESULTAT);
    float* V1s = (float*)malloc(MAX_SIZE*sizeof(float)); float* V2s = (float*)malloc(MAX_SIZE*sizeof(float)); float* V2s_copy = (float*)malloc(MAX_SIZE*sizeof(float));
    double* V1d = (double*)malloc(MAX_SIZE*sizeof(double)); double* V2d = (double*)malloc(MAX_SIZE*sizeof(double)); double* V2d_copy = (double*)malloc(MAX_SIZE*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)); complexe_float_t* V2c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)); complexe_float_t* V2c_copy = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)); complexe_double_t* V2z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)); complexe_double_t* V2z_copy = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t));
    
    float* Ms = (float*)malloc(M_RESULTAT*N_RESULTAT*sizeof(float));
    double* Md = (double*)malloc(M_RESULTAT*N_RESULTAT*sizeof(double));
    complexe_float_t* Mc = (complexe_float_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* Mz = (complexe_double_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_double_t));
    
    
    
    complexe_float_t tmpc = gen_complexe_float(2,0);
    complexe_double_t tmpz = gen_complexe_double(2,0);

    


    printf("                   TEST_GEMV\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("-------------------------- test 1 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor): \n");
    void_matrix_sinit(Ms,1,M_RESULTAT,N_RESULTAT);
    void_vector_sinit(V1s, 1,N_RESULTAT);
    void_vector_sinit(V2s, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 2 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor): \n");
    void_vector_sinit(V2s, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);


    printf("-------------------------- test 3 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor): \n");
    void_matrix_sinit_rand(Ms,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,N_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,M_RESULTAT);
    mncblas_scopy(M_RESULTAT,V2s,1,V2s_copy,1);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 4 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_FLOAT,Ms,M_RESULTAT,N_RESULTAT);
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s_copy,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,M_RESULTAT);


    printf("-------------------------- test 5 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_sinit(Ms,1,M_RESULTAT,N_RESULTAT);
    void_vector_sinit(V1s, 1,M_RESULTAT);
    void_vector_sinit(V2s, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 6 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_sinit(V2s, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);


    printf("-------------------------- test 7 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_sinit_rand(Ms,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,M_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,N_RESULTAT);
    mncblas_scopy(N_RESULTAT,V2s,1,V2s_copy,1);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 8 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_FLOAT,Ms,M_RESULTAT,N_RESULTAT);
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s_copy,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,N_RESULTAT);

    printf("<------------------------------------------>\n                   double\n");

    printf("-------------------------- test 1 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor): \n");
    void_matrix_dinit(Md,1,M_RESULTAT,N_RESULTAT);
    void_vector_dinit(V1d, 1,N_RESULTAT);
    void_vector_dinit(V2d, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 2 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor): \n");
    void_vector_dinit(V2d, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);


    printf("-------------------------- test 3 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor): \n");
    void_matrix_dinit_rand(Md,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,N_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,M_RESULTAT);
    mncblas_dcopy(M_RESULTAT,V2d,1,V2d_copy,1);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 4 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_DOUBLE,Md,M_RESULTAT,N_RESULTAT);
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d_copy,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,M_RESULTAT);


    printf("-------------------------- test 5 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_dinit(Md,1,M_RESULTAT,N_RESULTAT);
    void_vector_dinit(V1d, 1,M_RESULTAT);
    void_vector_dinit(V2d, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 6 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_dinit(V2d, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);


    printf("-------------------------- test 7 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_dinit_rand(Md,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,M_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,N_RESULTAT);
    mncblas_dcopy(N_RESULTAT,V2d,1,V2d_copy,1);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 8 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_DOUBLE,Md,M_RESULTAT,N_RESULTAT);
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d_copy,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,N_RESULTAT);



    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("-------------------------- test 1 (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor): \n");
    void_matrix_cinit(Mc,gen_complexe_float(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_cinit(V1c, gen_complexe_float(1,0),N_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(2,0),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 1_bis (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor): \n");
    void_matrix_cinit(Mc,gen_complexe_float(0,1),M_RESULTAT,N_RESULTAT);
    void_vector_cinit(V1c, gen_complexe_float(0,1),N_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(0,2),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 1_ter (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor): \n");
    void_matrix_cinit(Mc,gen_complexe_float(1,1),M_RESULTAT,N_RESULTAT);
    void_vector_cinit(V1c, gen_complexe_float(1,1),N_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(2,2),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 2 (V2c=2*Mc*V1c+2*V2c) (MNCblasColMajor): \n");
    void_matrix_cinit(Mc,gen_complexe_float(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_cinit(V1c, gen_complexe_float(1,0),N_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(2,0),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);


    printf("-------------------------- test 3 (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor): \n");
    void_matrix_cinit_rand(Mc,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,N_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,M_RESULTAT);
    mncblas_ccopy(M_RESULTAT,V2c,1,V2c_copy,1);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 4 (V2c=2*Mc*V1c+2*V2c) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_COMPLEXE_FLOAT,Mc,M_RESULTAT,N_RESULTAT);
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("V2c : "); vector_print(V2c_copy,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    mncblas_cgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c_copy,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c_copy,TYPE_COMPLEXE_FLOAT,M_RESULTAT);


    printf("-------------------------- test 5 (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_cinit(Mc,gen_complexe_float(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_cinit(V1c, gen_complexe_float(1,0),M_RESULTAT);
    void_vector_cinit(V2c, gen_complexe_float(2,0),N_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 6 (V2c=2*Mc*V1c+2*V2c) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_cinit(V2c, gen_complexe_float(2,0),N_RESULTAT);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    mncblas_cgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);


    printf("-------------------------- test 7 (V2c=2*Mc*V1c+2*V2c) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_cinit_rand(Mc,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,M_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,N_RESULTAT);
    mncblas_ccopy(N_RESULTAT,V2c,1,V2c_copy,1);
    printf("---- Avant :\n");
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    mncblas_cgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 8 (V2c=2*Mc*V1c+2*V2c) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_COMPLEXE_FLOAT,Mc,M_RESULTAT,N_RESULTAT);
    printf("Mc : "); matrix_print(Mc,TYPE_COMPLEXE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,M_RESULTAT);
    printf("V2c : "); vector_print(V2c_copy,TYPE_COMPLEXE_FLOAT,N_RESULTAT);
    mncblas_cgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpc,Mc,1,V1c,1,&tmpc,V2c_copy,1);
    printf("---- Apres :\n");
    printf("V2c : "); vector_print(V2c_copy,TYPE_COMPLEXE_FLOAT,N_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("-------------------------- test 1 (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor): \n");
    void_matrix_zinit(Mz,gen_complexe_double(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_zinit(V1z, gen_complexe_double(1,0),N_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(2,0),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 1_bis (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor): \n");
    void_matrix_zinit(Mz,gen_complexe_double(0,1),M_RESULTAT,N_RESULTAT);
    void_vector_zinit(V1z, gen_complexe_double(0,1),N_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(0,2),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 1_ter (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor): \n");
    void_matrix_zinit(Mz,gen_complexe_double(1,1),M_RESULTAT,N_RESULTAT);
    void_vector_zinit(V1z, gen_complexe_double(1,1),N_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(2,2),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 2 (V2z=2*Mz*V1z+2*V2z) (MNCblasColMajor): \n");
    void_matrix_zinit(Mz,gen_complexe_double(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_zinit(V1z, gen_complexe_double(1,0),N_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(2,0),M_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);


    printf("-------------------------- test 3 (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor): \n");
    void_matrix_zinit_rand(Mz,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,N_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,M_RESULTAT);
    mncblas_zcopy(M_RESULTAT,V2z,1,V2z_copy,1);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 4 (V2z=2*Mz*V1z+2*V2z) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_COMPLEXE_DOUBLE,Mz,M_RESULTAT,N_RESULTAT);
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("V2z : "); vector_print(V2z_copy,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    mncblas_zgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z_copy,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z_copy,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);


    printf("-------------------------- test 5 (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_zinit(Mz,gen_complexe_double(1,0),M_RESULTAT,N_RESULTAT);
    void_vector_zinit(V1z, gen_complexe_double(1,0),M_RESULTAT);
    void_vector_zinit(V2z, gen_complexe_double(2,0),N_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 6 (V2z=2*Mz*V1z+2*V2z) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_zinit(V2z, gen_complexe_double(2,0),N_RESULTAT);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    mncblas_zgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);


    printf("-------------------------- test 7 (V2z=2*Mz*V1z+2*V2z) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_zinit_rand(Mz,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,M_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,N_RESULTAT);
    mncblas_zcopy(N_RESULTAT,V2z,1,V2z_copy,1);
    printf("---- Avant :\n");
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    mncblas_zgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 8 (V2z=2*Mz*V1z+2*V2z) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_COMPLEXE_DOUBLE,Mz,M_RESULTAT,N_RESULTAT);
    printf("Mz : "); matrix_print(Mz,TYPE_COMPLEXE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,M_RESULTAT);
    printf("V2z : "); vector_print(V2z_copy,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);
    mncblas_zgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,&tmpz,Mz,1,V1z,1,&tmpz,V2z_copy,1);
    printf("---- Apres :\n");
    printf("V2z : "); vector_print(V2z_copy,TYPE_COMPLEXE_DOUBLE,N_RESULTAT);

    free(V1s);
    free(V1d);
    free(V1c);
    free(V1z);
    free(V2s);
    free(V2s_copy);
    free(V2d);
    free(V2d_copy);
    free(V2c);
    free(V2c_copy);
    free(V2z);
    free(V2z_copy);
    free(Ms);
    free(Md);
    free(Mc);
    free(Mz);

    #define M_FLOPS             1000
    #define N_FLOPS             1000
    #define NB_EXPE_VISIBLE     6
    #define NB_EXPE             100
    #define NB_OPE_REEL        (M_FLOPS*N_FLOPS*2+M_FLOPS*3)
    #define NB_OPE_COMPLEXE    (M_FLOPS*N_FLOPS*8+M_FLOPS*14)

    MAX_SIZE = (M_RESULTAT < N_RESULTAT ? N_RESULTAT : M_RESULTAT);
    V1s = (float*)malloc(MAX_SIZE*sizeof(float)); V2s = (float*)malloc(MAX_SIZE*sizeof(float));
    V1d = (double*)malloc(MAX_SIZE*sizeof(double)); V2d = (double*)malloc(MAX_SIZE*sizeof(double));
    V1c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)); V2c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)); V2z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t));

    Ms = (float*)malloc(M_FLOPS*N_FLOPS*sizeof(float));
    Md = (double*)malloc(M_FLOPS*N_FLOPS*sizeof(double));
    Mc = (complexe_float_t*)malloc(M_FLOPS*N_FLOPS*sizeof(complexe_float_t));
    Mz = (complexe_double_t*)malloc(M_FLOPS*N_FLOPS*sizeof(complexe_double_t));

    init_flop();

    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Ms,1,V1s,1,2,V2s,1);
        end = _rdtsc();
        calcul_flop("mncblas_sgemv : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Ms,1,V1s,1,2,V2s,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_sgemv : ", NB_EXPE*NB_OPE_REEL ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Md,1,V1d,1,2,V2d,1);
        end = _rdtsc();
        calcul_flop("mncblas_dgemv : ", NB_OPE_REEL ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Md,1,V1d,1,2,V2d,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_dgemv : ", NB_EXPE*(NB_OPE_REEL) ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
        end = _rdtsc();
        calcul_flop("mncblas_cgemv : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,&tmpc,Mc,1,V1c,1,&tmpc,V2c,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cgemv : ", NB_EXPE*(NB_OPE_COMPLEXE) ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
        end = _rdtsc();
        calcul_flop("mncblas_zgemv : ", NB_OPE_COMPLEXE ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,&tmpz,Mz,1,V1z,1,&tmpz,V2z,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zgemv : ", NB_EXPE*(NB_OPE_COMPLEXE) ,end-start);


    // free(V1s);
    // free(V1d);
    // free(V1c);
    // free(V1z);
    // free(V2s);
    // free(V2d);
    // free(V2c);
    // free(V2z);
    // free(Ms);
    // free(Md);
    // free(Mc);
    // free(Mz);
}