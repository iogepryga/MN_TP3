#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

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
  for (register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    dot += X [i] * Y [j];
  }
  return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY) { // NB OPE FLOTANTE = 2*N
  register double dot = 0.0;
  for (register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    dot += X [i] * Y [j];
  }
  return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu) { // NB OPE FLOTANTE = 8*N
  register complexe_float_t dot; dot.real=0; dot.imaginary =0; //tmp
  for(register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    //tentative 1 :
    /*tmp = mult_complexe_float(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+j));
    dot.real += tmp.real;
    dot.imaginary += tmp.imaginary;*/

    //tentative 2 :
    dot = add_complexe_float(dot,mult_complexe_float(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+j))); // 6 + 2
  }
  if(N >= 0) {
    *(complexe_float_t*)dotu = dot;
  } else {
    ((complexe_float_t*)dotu)->real = 0;((complexe_float_t*)dotu)->imaginary = 0;
  }
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc) { // NB OPE FLOTANTE = 9*N
  register complexe_float_t dot; dot.real=0; dot.imaginary =0; //tmp
  for(register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    //tentative 1 :
    /*tmp = mult_complexe_double(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+j));
    dot.real += tmp.real;
    dot.imaginary += tmp.imaginary;*/

    //tentative 2 :
    ((complexe_float_t*)X+i)->imaginary *= -1;
    dot = add_complexe_float(dot,mult_complexe_float(*((complexe_float_t*)X+i),*((complexe_float_t*)Y+j)));
  }
  if(N >= 0) {
    *(complexe_float_t*)dotc = dot;
  } else {
    ((complexe_float_t*)dotc)->real = 0;((complexe_float_t*)dotc)->imaginary = 0;
  }
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu) { // NB OPE FLOTANTE = 8*N
  register complexe_double_t dot; dot.real=0; dot.imaginary =0; //tmp
  for(register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    dot = add_complexe_double(dot,mult_complexe_double(*((complexe_double_t*)X+i),*((complexe_double_t*)Y+j)));
  }
  if(N >= 0) {
    *(complexe_double_t*)dotu = dot;
  } else {
    ((complexe_double_t*)dotu)->real = 0;((complexe_double_t*)dotu)->imaginary = 0;
  }
}
  
void mncblas_zdotc_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc) { // NB OPE FLOTANTE = 9*N
  register complexe_double_t dot; dot.real=0; dot.imaginary =0; //tmp
  for(register unsigned int i = 0 , j= 0 ; i < N ; i += incX , j+=incY) {
    ((complexe_double_t*)X+i)->imaginary *= -1;
    dot = add_complexe_double(dot,mult_complexe_double(*((complexe_double_t*)X+i),*((complexe_double_t*)Y+j)));
  }
  if(N >= 0) {
    *(complexe_double_t*)dotc = dot;
  } else {
    ((complexe_double_t*)dotc)->real = 0;((complexe_double_t*)dotc)->imaginary = 0;
  }
}




