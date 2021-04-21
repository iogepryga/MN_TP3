#include "mnblas.h"
#include "complexe.h"

CBLAS_INDEX mnblas_isamin(const int N, const float *X, const int incX) {
    register CBLAS_INDEX imin;
    register float min;
    if(N <= 0) {
        return 0;
    } else {
        imin = 0;
        min = X[0];
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        if(X[i] < min) {
            min = X[i];
            imin = i;
        }
    }
    return imin;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX) {
    register CBLAS_INDEX imin;
    register double min;
    if(N <= 0) {
        return 0;
    } else {
        imin = 0;
        min = X[0];
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        if(X[i] < min) {
            min = X[i];
            imin = i;
        }
    }
    return imin;
}

CBLAS_INDEX mnblas_icamin(const int N, const void *X, const int incX) {
    register CBLAS_INDEX imin;
    register float min,tmp;
    if(N <= 0) {
        return 0;
    } else {
        imin = 0;
        min = (((complexe_float_t*)X)->real < 0 ? -((complexe_float_t*)X)->real : ((complexe_float_t*)X)->real)
        + (((complexe_float_t*)X)->imaginary < 0 ? -((complexe_float_t*)X)->imaginary : ((complexe_float_t*)X)->imaginary);
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        tmp = (((complexe_float_t*)X+i)->real < 0 ? -((complexe_float_t*)X+i)->real : ((complexe_float_t*)X+i)->real)
        + (((complexe_float_t*)X+i)->imaginary < 0 ? -((complexe_float_t*)X+i)->imaginary : ((complexe_float_t*)X+i)->imaginary);
        if(tmp < min) {
            min = tmp;
            imin = i;
        }
    }
    return imin;
}

CBLAS_INDEX mnblas_izamin(const int N, const void *X, const int incX) {
    register CBLAS_INDEX imin;
    register double min,tmp;
    if(N <= 0) {
        return 0;
    } else {
        imin = 0;
        min = (((complexe_double_t*)X)->real < 0 ? -((complexe_double_t*)X)->real : ((complexe_double_t*)X)->real)
        + (((complexe_double_t*)X)->imaginary < 0 ? -((complexe_double_t*)X)->imaginary : ((complexe_double_t*)X)->imaginary);
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        tmp = (((complexe_double_t*)X+i)->real < 0 ? -((complexe_double_t*)X+i)->real : ((complexe_double_t*)X+i)->real)
        + (((complexe_double_t*)X+i)->imaginary < 0 ? -((complexe_double_t*)X+i)->imaginary : ((complexe_double_t*)X+i)->imaginary);
        if(tmp < min) {
            min = tmp;
            imin = i;
        }
    }
    return imin;
}