#include "mnblas.h"
#include "complexe.h"


void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const float alpha, const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY) { // NB OPE FLOTANTE = 2*M*N + 3*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register float sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum = 0;
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+i*N+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sum ;
            for(register unsigned int i = 0; i < N ; i+= incY) {
                sum = 0;
                for(register int j = 0; j < M ; j+= incX) {
                    sum += *(A+j*N+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register float sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum = 0;
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+j*M+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum = 0;
                for(j = 0; j < M ; j+= incX) {
                    sum += *(A+i*M+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY) { // NB OPE FLOTANTE = 2*M*N + 3*M
   if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register double sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum = 0;
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+i*N+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register double sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum = 0;
                for(j = 0; j < M ; j+= incX) {
                    sum += *(A+j*N+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register double sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum = 0;
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+j*M+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register double sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum = 0;
                for(j = 0; j < M ; j+= incX) {
                    sum += *(A+i*M+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) { // NB OPE FLOTANTE = 8*M*N + 14*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+i*N+j)*X[j];
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*N+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+j*N+i)*X[j];
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+j*N+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+j*N+i)*X[j];
                    (((complexe_float_t*)A)+j*N+i)->imaginary *= -1;
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+j*N+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+j*M+i)*X[j];
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+j*M+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+i*M+j)*X[j];
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*M+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
        else if (TransA == MNCblasConjTrans) {
            register complexe_float_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+i*M+j)*X[j];
                    (((complexe_float_t*)A)+i*M+j)->imaginary *= -1;
                    sum = add_complexe_float(sum,mult_complexe_float(*(((complexe_float_t*)A)+i*M+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sum),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY){ // NB OPE FLOTANTE = 10*M*N + 14*M
    /*if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sum ;
            for(register int i = 0, j; i < M ; i+= incY) {
                 sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+i*N+j)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sum ;
            for(register int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+j*N+i)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+j*M+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sum ;
            for(register int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+j*M+i)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+j*M+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sum ;
            for(register int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+i*M+j)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    }*/

    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+i*N+j)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*N+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+j*N+i)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+j*N+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+j*N+i)*X[j];
                    (((complexe_double_t*)A)+j*N+i)->imaginary *= -1;
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+j*N+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    //sum += *(A+j*M+i)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+j*M+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+i*M+j)*X[j];
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*M+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_double_t sum ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sum.real = 0; sum.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    //sum += *(A+i*M+j)*X[j];
                    (((complexe_double_t*)A)+i*M+j)->imaginary *= -1;
                    sum = add_complexe_double(sum,mult_complexe_double(*(((complexe_double_t*)A)+i*M+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sum),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    }
}