#include "mnblas.h"
#include "complexe.h"
#include "math.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY) { // NB OPE FLOTANTE = 2*N
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0 ; i < max ; i++) {
#pragma omp critical
    Y[diff*i] += alpha * X[i];
  }
}

void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY) { // NB OPE FLOTANTE = 2*N
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0 ; i < max ; i++) {
#pragma omp critical
    Y[diff*i] += alpha * X[i];
	}
}

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) { // NB OPE FLOTANTE = (6 + 2) * N = 8*N
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0 ; i < max ; i++) {
#pragma omp critical
    *((complexe_float_t*)Y+diff*i) = add_complexe_float(*((complexe_float_t*)Y+diff*i),mult_complexe_float(*((complexe_float_t*)alpha),*((complexe_float_t*)X+i)));
  }
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) { // NB OPE FLOTANTE = (6 + 2) * N = 8*N
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0 ; i < max ; i++) {
#pragma omp critical
    *((complexe_double_t*)Y+diff*i) = add_complexe_double(*((complexe_double_t*)Y+diff*i),mult_complexe_double(*((complexe_double_t*)alpha),*((complexe_double_t*)X+i)));
  }
}