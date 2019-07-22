// calculates the first partial derivative of the discretized action to second order  with respect to the ij-th variable

#ifndef partderij2h
#define partderij2h


void parder2ij(mpfr_t** a, mpfr_t * par, int N, int m, int o, int i, int j, mpfr_t J,mpfr_t** M, mpfr_t zero, int prec);

#endif
