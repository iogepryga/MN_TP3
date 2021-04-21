#include "mnblas.h"
#include "complexe.h"

CBLAS_INDEX mnblas_isamax(const int N, const float *X, const int incX) {
    register CBLAS_INDEX imax;
    register float max;
    if(N <= 0) {
        return 0;
    } else {
        imax = 0;
        max = X[0];
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        if(X[i] > max) {
            max = X[i];
            imax = i;
        }
    }
    return imax;
}

CBLAS_INDEX mnblas_idamax(const int N, const double *X, const int incX) {
    register CBLAS_INDEX imax;
    register double max;
    if(N <= 0) {
        return 0;
    } else {
        imax = 0;
        max = X[0];
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        if(X[i] > max) {
            max = X[i];
            imax = i;
        }
    }
    return imax;
}

CBLAS_INDEX mnblas_icamax(const int N, const void *X, const int incX) {
    register CBLAS_INDEX imax;
    register float max,tmp;
    if(N <= 0) {
        return 0;
    } else {
        imax = 0;
        max = (((complexe_float_t*)X)->real < 0 ? -((complexe_float_t*)X)->real : ((complexe_float_t*)X)->real)
        + (((complexe_float_t*)X)->imaginary < 0 ? -((complexe_float_t*)X)->imaginary : ((complexe_float_t*)X)->imaginary);
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        tmp = (((complexe_float_t*)X+i)->real < 0 ? -((complexe_float_t*)X+i)->real : ((complexe_float_t*)X+i)->real)
        + (((complexe_float_t*)X+i)->imaginary < 0 ? -((complexe_float_t*)X+i)->imaginary : ((complexe_float_t*)X+i)->imaginary);
        if(tmp > max) {
            max = tmp;
            imax = i;
        }
    }
    return imax;
}

CBLAS_INDEX mnblas_izamax(const int N, const void *X, const int incX) {
    register CBLAS_INDEX imax;
    register double max,tmp;
    if(N <= 0) {
        return 0;
    } else {
        imax = 0;
        max = (((complexe_double_t*)X)->real < 0 ? -((complexe_double_t*)X)->real : ((complexe_double_t*)X)->real)
        + (((complexe_double_t*)X)->imaginary < 0 ? -((complexe_double_t*)X)->imaginary : ((complexe_double_t*)X)->imaginary);
    }
    for (register unsigned int i = 1; i < N ; i += incX) {
        tmp = (((complexe_double_t*)X+i)->real < 0 ? -((complexe_double_t*)X+i)->real : ((complexe_double_t*)X+i)->real)
        + (((complexe_double_t*)X+i)->imaginary < 0 ? -((complexe_double_t*)X+i)->imaginary : ((complexe_double_t*)X+i)->imaginary);
        if(tmp > max) {
            max = tmp;
            imax = i;
        }
    }
    return imax;
}