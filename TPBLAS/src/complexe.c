#include "complexe.h"

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2) {
	complexe_float_t r;
	r.real = c1.real + c2.real;
	r.imaginary = c1.imaginary + c2.imaginary;
	return r;
}

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2) {
	complexe_double_t r;
	r.real = c1.real + c2.real;
	r.imaginary = c1.imaginary + c2.imaginary;
	return r;
}

complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;
  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
  r.imaginary = c1.imaginary * c2.real + c1.real * c2.imaginary;
  return r ;
}

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
  {
  complexe_double_t r ;
  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
  r.imaginary = c1.imaginary * c2.real + c1.real * c2.imaginary;
  return r ;
}

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2) {
  	/*
	c1 = a + ib
  	c2 = c + id
  	c1 / c2 =   (c1 * c2_ ) / (c2 * c2_)
		  =   (c1 * c2_ ) / (c*c + d*d)
		  =   (a + ib) * (c - id) / (c² + d²)
	*/
	complexe_float_t r ;
	r.real = c2.real;
	r.imaginary = -c2.imaginary;
	float br = c2.real*c2.real + c1.imaginary*c1.imaginary;
	r = mult_complexe_float(c1,r);
	r.real /= br;
  	r.imaginary /= br;
	return r;
}

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2) {
	complexe_double_t r;
	r.real = c2.real;
	r.imaginary = -c2.imaginary;
	double br = c2.real*c2.real + c1.imaginary*c1.imaginary;
	r = mult_complexe_double(c1,r);
	r.real /= br;
  	r.imaginary /= br;
	return r;
}
	
