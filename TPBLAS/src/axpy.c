#include "mnblas.h"
#include "complexe.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY) { // NB OPE FLOTANTE = 2*N
    for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
        Y[j] += alpha * X[i];
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY) { // NB OPE FLOTANTE = 2*N
    for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
        Y[j] += alpha * X[i];
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) { // NB OPE FLOTANTE = (6 + 2) * N = 8*N
    for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
        *((complexe_float_t*)Y+j) = add_complexe_float(*((complexe_float_t*)Y+j),mult_complexe_float(*((complexe_float_t*)alpha),*((complexe_float_t*)X+i)));
    }
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) { // NB OPE FLOTANTE = (6 + 2) * N = 8*N
    for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
        *((complexe_double_t*)Y+j) = add_complexe_double(*((complexe_double_t*)Y+j),mult_complexe_double(*((complexe_double_t*)alpha),*((complexe_double_t*)X+i)));
    }
}