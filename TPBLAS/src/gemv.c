#include "mnblas.h"
#include "complexe.h"


void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const float alpha, const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY) { // NB OPE FLOTANTE = 2*M*N + 3*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans) {
#pragma omp parallel for
        for(register unsigned int i = 0; i < M ; i+= incY) {
            register float sum = 0;
            for(register unsigned int j = 0; j < N ; j+= incX) {
                sum += *(A+i*N+j)*X[j];
            }
            Y[i] = alpha*sum + beta*Y[i];
        }
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY) { // NB OPE FLOTANTE = 2*M*N + 3*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans) {
#pragma omp parallel for
        for(register unsigned int i = 0; i < M ; i+= incY) {
            register double sum = 0;
            for(register unsigned int j = 0; j < N ; j+= incX) {
                sum += *(A+i*N+j)*X[j];
            }
            Y[i] = alpha*sum + beta*Y[i];
        }
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) { // NB OPE FLOTANTE = 8*M*N + 14*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans) {
#pragma omp parallel for
        for(register unsigned int i = 0; i < M ; i+= incY) {
            register complexe_float_t sum ; sum.real = 0; sum.imaginary = 0;
            for(register unsigned int j = 0; j < N ; j+= incX) {
                //sum += *(A+i*N+j)*X[j];
                sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*N+j),*(((complexe_float_t*)X)+j)));
            }
            *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY){ // NB OPE FLOTANTE = 10*M*N + 14*M
    if(layout == MNCblasRowMajor && TransA == MNCblasNoTrans) {
#pragma omp parallel for
        for(register unsigned int i = 0; i < M ; i+= incY) {
            register complexe_double_t sum; sum.real = 0; sum.imaginary = 0;
            for(register unsigned int j = 0; j < N ; j+= incX) {
                //sum += *(A+i*N+j)*X[j];
                sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+j),*(((complexe_double_t*)X)+j)));
            }
            *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
        }
    } 
}