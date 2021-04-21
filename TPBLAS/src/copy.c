#include "mnblas.h"
#include "complexe.h"

void mncblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY) { // N AFFECTIONS -> N*sizeof(float) octets
  for (register unsigned int i = 0, j = 0 ; ((i < N) && (j < N)) ; i += incX, j += incY){
    Y [j] = X [i] ;
  }
}

void mncblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY) { // N AFFECTIONS -> N*sizeof(double) octets
  for (register unsigned int i = 0 , j = 0 ; ((i < N) && (j < N)) ; i += incX, j += incY){
    Y [j] = X [i] ;
  }
}

void mncblas_ccopy(const int N, const void *X, const int incX, void *Y, const int incY) { // N AFFECTIONS -> N*sizeof(complexe_float_t) octets
  for (register unsigned int i = 0 , j = 0 ; ((i < N) && (j < N)) ; i += incX, j += incY){
    *((complexe_float_t*)Y+j) = *((complexe_float_t*)X+i) ;
  }
}

void mncblas_zcopy(const int N, const void *X, const int incX, void *Y, const int incY) { // N AFFECTIONS -> N*sizeof(complexe_double_t) octets
  for (register unsigned int i = 0 , j = 0 ; ((i < N) && (j < N)) ; i += incX, j += incY){
    *((complexe_double_t*)Y+j) = *((complexe_double_t*)X+i) ;
  }
}

