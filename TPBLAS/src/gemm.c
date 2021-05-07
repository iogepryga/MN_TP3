#include "mnblas.h"
#include "complexe.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc) { // NB OPE FLOTANTE = 2*M*N*K + 3*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans && TransB == MNCblasNoTrans) {
#pragma omp parallel for
        for(register int i = 0; i < M ; i++) {
            for(register int j = 0; j < K ; j++) {
                register float sum = 0;
                for(register int k = 0; k < N ;k++)
                    sum += *(A+i*N+k) * B[k*K+j];
                C[i*K+j] = alpha*sum + beta*C[i*K+j];
             }
        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) { // NB OPE FLOTANTE = 2*M*N*K + 3*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans && TransB == MNCblasNoTrans) {
#pragma omp parallel for
        for(register int i = 0; i < M ; i++) {
            for(register int j = 0; j < K ; j++) {
                register double sum = 0;
                for(register int k = 0; k < N ;k++)
                    sum += *(A+i*N+k) * B[k*K+j];
                C[i*K+j] = alpha*sum + beta*C[i*K+j];
            }
         }
    }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc) { // NB OPE FLOTANTE = 10*M*N*K + 18*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans && TransB == MNCblasNoTrans) {
#pragma omp parallel for
        for(register int i = 0; i < M ; i++) {
            for(register int j = 0; j < K ; j++) {
                register complexe_float_t sum ; sum.real = 0; sum.imaginary = 0;
                for(register int k = 0; k < N ;k++) {
                    sum.real += (((complexe_float_t*)A)+i*N+k)->real * (((complexe_float_t*)B)+k*K+j)->real - (((complexe_float_t*)A)+i*N+k)->imaginary * (((complexe_float_t*)B)+k*K+j)->imaginary;
                    sum.imaginary += (((complexe_float_t*)A)+i*N+k)->imaginary * (((complexe_float_t*)B)+k*K+j)->real + (((complexe_float_t*)A)+i*N+k)->real * (((complexe_float_t*)B)+k*K+j)->imaginary;
                    // sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*N+k),*(((complexe_float_t*)B)+k*K+j)));
                }
                (((complexe_float_t*)C)+i*K+j)->real = ((complexe_float_t*)alpha)->real * sum.real - ((complexe_float_t*)alpha)->imaginary * sum.imaginary + ((complexe_float_t*)beta)->real * (((complexe_float_t*)C)+i*K+j)->real - ((complexe_float_t*)beta)->imaginary * (((complexe_float_t*)C)+i*K+j)->imaginary;
                (((complexe_float_t*)C)+i*K+j)->imaginary = ((complexe_float_t*)alpha)->imaginary * sum.real + ((complexe_float_t*)alpha)->real * sum.imaginary + ((complexe_float_t*)beta)->imaginary * (((complexe_float_t*)C)+i*K+j)->real + ((complexe_float_t*)beta)->real * (((complexe_float_t*)C)+i*K+j)->imaginary;
                // *(((complexe_float_t*)C)+i*K+j) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)C)+i*K+j)));
            }
        }
    }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){ // NB OPE FLOTANTE = 10*M*N*K + 18*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans && TransB == MNCblasNoTrans) {
#pragma omp parallel for
        for(register int i = 0; i < M ; i++) {
            for(register int j = 0; j < K ; j++) {
                register complexe_double_t sum ; sum.real = 0; sum.imaginary = 0;
                for(register int k = 0; k < N ;k++) {
                    sum.real += (((complexe_double_t*)A)+i*N+k)->real * (((complexe_double_t*)B)+k*K+j)->real - (((complexe_double_t*)A)+i*N+k)->imaginary * (((complexe_double_t*)B)+k*K+j)->imaginary;
                    sum.imaginary += (((complexe_double_t*)A)+i*N+k)->imaginary * (((complexe_double_t*)B)+k*K+j)->real + (((complexe_double_t*)A)+i*N+k)->real * (((complexe_double_t*)B)+k*K+j)->imaginary;
                    // sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+k),*(((complexe_double_t*)B)+k*K+j)));
                }
                (((complexe_float_t*)C)+i*K+j)->real = ((complexe_double_t*)alpha)->real * sum.real - ((complexe_double_t*)alpha)->imaginary * sum.imaginary + ((complexe_double_t*)beta)->real * (((complexe_double_t*)C)+i*K+j)->real - ((complexe_double_t*)beta)->imaginary * (((complexe_double_t*)C)+i*K+j)->imaginary;
                (((complexe_float_t*)C)+i*K+j)->imaginary = ((complexe_double_t*)alpha)->imaginary * sum.real + ((complexe_double_t*)alpha)->real * sum.imaginary + ((complexe_double_t*)beta)->imaginary * (((complexe_double_t*)C)+i*K+j)->real + ((complexe_double_t*)beta)->real * (((complexe_double_t*)C)+i*K+j)->imaginary;
                // *(((complexe_double_t*)C)+i*K+j) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)C)+i*K+j)));
            }
        }
    }
}