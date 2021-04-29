#include "mnblas.h"
#include "complexe.h"
#include <math.h>

void mncblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY) { // N AFFECTIONS -> N*sizeof(float) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0; i < max ; i++) {
    Y [diff*i] = X [i];
  }
}

void mncblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY) { // N AFFECTIONS -> N*sizeof(double) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0; i < max ; i++) {
    Y [diff*i] = X [i] ;
  }
}

void mncblas_ccopy(const int N, const void *X, const int incX, void *Y, const int incY) { // N AFFECTIONS -> N*sizeof(complexe_float_t) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0; i < max ; i++) {
    *((complexe_float_t*)Y+diff*i) = *((complexe_float_t*)X+i) ;
  }
}

void mncblas_zcopy(const int N, const void *X, const int incX, void *Y, const int incY) { // N AFFECTIONS -> N*sizeof(complexe_double_t) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for
  for (register unsigned int i = 0; i < max ; i++) {
    *((complexe_double_t*)Y+diff*i) = *((complexe_double_t*)X+i) ;
  }
}

