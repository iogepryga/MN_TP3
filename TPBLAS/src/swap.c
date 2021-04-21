#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX, float *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(float) octets
  register float save;
  for (register unsigned int i = 0, j= 0; ((i < N) && (j < N)) ; i += incX, j+=incY) {
    save = Y [j];
    Y [j] = X [i];
    X [i] = save;
  }
}

void mncblas_dswap(const int N, double *X, const int incX, double *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(double) octets
  register double save;
  for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
    save = Y [j];
    Y [j] = X [i];
    X [i] = save;
  }
}

void mncblas_cswap(const int N, void *X, const int incX, void *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(complexe_float_t) octets
  register complexe_float_t save ;
  for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
    save = *((complexe_float_t*)Y+j);
    *((complexe_float_t*)Y+j) = *((complexe_float_t*)X+i);
    *((complexe_float_t*)X+i) = save;
  }
}

void mncblas_zswap(const int N, void *X, const int incX, void *Y, const int incY) { // 3*N AFFECTIONS -> 3*N*sizeof(complexe_double_t) octets
  register complexe_double_t save;
  for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j+=incY) {
    save = *((complexe_double_t*)Y+j);
    *((complexe_double_t*)Y+j) = *((complexe_double_t*)X+i);
    *((complexe_double_t*)X+i) = save;
  }
}

