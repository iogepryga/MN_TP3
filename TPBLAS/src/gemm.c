#include "mnblas.h"
#include "complexe.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc) { // NB OPE FLOTANTE = 2*M*N*K + 3*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register float sum ;
            for(register int i = 0; i < M ; i++) {
                for(register int j = 0; j < K ; j++) {
                    sum = 0;
                    for(register int k = 0; k < N ;k++) {
                        sum += *(A+i*N+k) * B[k*K+j];
                    }
                    C[i*K+j] = alpha*sum + beta*C[i*K+j];
                }
            }
        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    } else if (layout == MNCblasColMajor) {
        if(TransA == MNCblasNoTrans) {

        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) { // NB OPE FLOTANTE = 2*M*N*K + 3*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register double sum ;
            for(register int i = 0; i < M ; i++) {
                for(register int j = 0; j < K ; j++) {
                    sum = 0;
                    for(register int k = 0; k < N ;k++) {
                        sum += *(A+i*N+k) * B[k*K+j];
                    }
                    C[i*K+j] = alpha*sum + beta*C[i*K+j];
                }
            }
        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    } else if (layout == MNCblasColMajor) {
        if(TransA == MNCblasNoTrans) {

        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc) { // NB OPE FLOTANTE = 10*M*N*K + 18*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_float_t sum ;
            for(register int i = 0; i < M ; i++) {
                for(register int j = 0; j < K ; j++) {
                    sum.real = 0; sum.imaginary = 0;
                    for(register int k = 0; k < N ;k++) {
                        //sum += *(A+i*N+k) * B[k*K+j];
                        sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*N+k),*(((complexe_float_t*)B)+k*K+j)));
                    }
                    //C[i*K+j] = alpha*sum + beta*C[i*K+j];
                    *(((complexe_float_t*)C)+i*K+j) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)C)+i*K+j)));
                }
            }
        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    } else if (layout == MNCblasColMajor) {
        if(TransA == MNCblasNoTrans) {

        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB,
                     const int M, const int N, const int K,
                     const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){ // NB OPE FLOTANTE = 10*M*N*K + 18*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sum ;
            for(register int i = 0; i < M ; i++) {
                for(register int j = 0; j < K ; j++) {
                    sum.real = 0; sum.imaginary = 0;
                    for(register int k = 0; k < N ;k++) {
                        //sum += *(A+i*N+k) * B[k*K+j];
                        sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+k),*(((complexe_double_t*)B)+k*K+j)));
                    }
                    //C[i*K+j] = alpha*sum + beta*C[i*K+j];
                    *(((complexe_double_t*)C)+i*K+j) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)C)+i*K+j)));
                }
            }
        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    } else if (layout == MNCblasColMajor) {
        if(TransA == MNCblasNoTrans) {

        } else if (TransA == MNCblasTrans) {

        } else if (TransA == MNCblasConjTrans) {

        }
    }
}