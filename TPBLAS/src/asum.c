#include "mnblas.h"
#include "complexe.h"

float  mnblas_sasum(const int N, const float *X, const int incX) { // NB OPE FLOTANTE = (2 MAX 1 MIN)*N
    register float sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += X[i] < 0 ? -X[i] : X[i];
    }
    return sum;
}

double mnblas_dasum(const int N, const double *X, const int incX) { // NB OPE FLOTANTE = (2 MAX 1 MIN)*N
    register double sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += X[i] < 0 ? -X[i] : X[i];
    }
    return sum;
}

float  mnblas_scasum(const int N, const void *X, const int incX) { // NB OPE FLOTANTE = (4 MAX 2 MIN)*N
    register float sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += (((complexe_float_t*)X+i)->real < 0 ? -((complexe_float_t*)X+i)->real : ((complexe_float_t*)X+i)->real)
        + (((complexe_float_t*)X+i)->imaginary < 0 ? -((complexe_float_t*)X+i)->imaginary : ((complexe_float_t*)X+i)->imaginary);
    }
    return sum;
}

double mnblas_dzasum(const int N, const void *X, const int incX) { // NB OPE FLOTANTE = (4 MAX 2 MIN)*N
    register double sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += (((complexe_double_t*)X+i)->real < 0 ? -((complexe_double_t*)X+i)->real : ((complexe_double_t*)X+i)->real)
        + (((complexe_double_t*)X+i)->imaginary < 0 ? -((complexe_double_t*)X+i)->imaginary : ((complexe_double_t*)X+i)->imaginary);
    }
    return sum;
}