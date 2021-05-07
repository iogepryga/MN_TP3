#include "mnblas.h"
#include "complexe.h"
#include <math.h>

void mncblas_sswap(const int N, float *X, const int incX, float *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(float) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
  register float save;
#pragma omp parallel for private(save)
  for (register unsigned int i = 0; i < max ; i++) {
    save = Y [i*diff];
    Y [i*diff] = X [i];
    X [i] = save;
  }
}

void mncblas_dswap(const int N, double *X, const int incX, double *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(double) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
  register double save;
#pragma omp parallel for private(save)
  for (register unsigned int i = 0; i < max ; i++) {
    save = Y [i*diff];
    Y [i*diff] = X [i];
    X [i] = save;
  }
}

void mncblas_cswap(const int N, void *X, const int incX, void *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(complexe_float_t) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
  register complexe_float_t save;
#pragma omp parallel for private(save)
  for (register unsigned int i = 0; i < max ; i++) {
    save = *((complexe_float_t*)Y+i*diff);
    *((complexe_float_t*)Y+i*diff) = *((complexe_float_t*)X+i);
    *((complexe_float_t*)X+i) = save;
  }
}

void mncblas_zswap(const int N, void *X, const int incX, void *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(complexe_double_t) octets
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
  register complexe_double_t save;
#pragma omp parallel for private(save)
  for (register unsigned int i = 0; i < max ; i++) {
    save = *((complexe_double_t*)Y+i*diff);
    *((complexe_double_t*)Y+i*diff) = *((complexe_double_t*)X+i);
    *((complexe_double_t*)X+i) = save;
  }
}

