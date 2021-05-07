#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"
#include "math.h"

/*
float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX, const float *Y, const int incY) { // NB OPE FLOTANTE = 2*N
  register float dot = 0.0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for reduction(+:dot)
  for (register unsigned int i = 0 ; i < max ; i++)
    dot += X [i] * Y [diff*i];
  return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY) { // NB OPE FLOTANTE = 2*N
  register double dot = 0.0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for reduction(+:dot)
  for (register unsigned int i = 0 ; i < max ; i++)
    dot += X [i] * Y [diff*i];
  return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu) { // NB OPE FLOTANTE = 8*N
  register complexe_float_t tmp;
  register float dotreal=0, dotimaginary=0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for private(tmp) reduction(+:dotreal) reduction(+:dotimaginary)
  for (register unsigned int i = 0 ; i < max ; i++) {
    tmp = mult_complexe_float(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+diff*i));
    dotreal +=  tmp.real;
    dotimaginary +=  tmp.imaginary;
  }
  
  if(N >= 0) {
    ((complexe_float_t*)dotu)->real = dotreal;
    ((complexe_float_t*)dotu)->imaginary = dotimaginary;
  } else {
    ((complexe_float_t*)dotu)->real = 0;
    ((complexe_float_t*)dotu)->imaginary = 0;
  }
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc) { // NB OPE FLOTANTE = 9*N
  register complexe_float_t tmp;
  register float dotreal=0, dotimaginary =0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for private(tmp) reduction(+:dotreal) reduction(+:dotimaginary) // on peut pas reduire une structure
  for (register unsigned int i = 0 ; i < max ; i++) {
    ((complexe_float_t*)X+i)->imaginary *= -1;
    tmp = mult_complexe_float(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+diff*i));
    dotreal +=  tmp.real;
    dotimaginary +=  tmp.imaginary;
  }

  if(N >= 0) {
    ((complexe_float_t*)dotc)->real = dotreal;
    ((complexe_float_t*)dotc)->imaginary = dotimaginary;
  } else {
    ((complexe_float_t*)dotc)->real = 0;
    ((complexe_float_t*)dotc)->imaginary = 0;
  }
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu) { // NB OPE FLOTANTE = 8*N
  register complexe_double_t tmp;
  register double dotreal=0, dotimaginary =0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for private(tmp) reduction(+:dotreal) reduction(+:dotimaginary) // on peut pas reduire une structure
  for (register unsigned int i = 0 ; i < max ; i++) {
    tmp = mult_complexe_double(*((complexe_double_t*)X+i),*((complexe_double_t*)Y+diff*i));
    dotreal +=  tmp.real;
    dotimaginary +=  tmp.imaginary;
  }

  if(N >= 0) {
    ((complexe_double_t*)dotu)->real = dotreal;
    ((complexe_double_t*)dotu)->imaginary = dotimaginary;
  } else {
    ((complexe_double_t*)dotu)->real = 0;
    ((complexe_double_t*)dotu)->imaginary = 0;
  }
}
  
void mncblas_zdotc_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc) { // NB OPE FLOTANTE = 9*N
  register complexe_double_t tmp;
  register double dotreal=0, dotimaginary =0;
  const int max = (incX < incY) ? ceil((float)N/(float)incY) : ceil((float)N/(float)incX);
  const int diff = (incY - incX + 1);
#pragma omp parallel for private(tmp) reduction(+:dotreal) reduction(+:dotimaginary) // on peut pas reduire une structure
  for (register unsigned int i = 0 ; i < max ; i++) {
    ((complexe_double_t*)X+i)->imaginary *= -1;
    tmp = mult_complexe_double(*((complexe_double_t*)X+i),*((complexe_double_t*)Y+diff*i));
    dotreal +=  tmp.real;
    dotimaginary +=  tmp.imaginary;
  }

  if(N >= 0) {
    ((complexe_double_t*)dotc)->real = dotreal;
    ((complexe_double_t*)dotc)->imaginary = dotimaginary;
  } else {
    ((complexe_double_t*)dotc)->real = 0;
    ((complexe_double_t*)dotc)->imaginary = 0;
  }
}




