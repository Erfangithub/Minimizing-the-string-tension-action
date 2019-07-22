// Takes N (number of gauge group), k (N-ality), m (number of partitions)
// , o (order), J (boundary float), prec (precision of calculation)
//
// Gives out a two dimensional array:
// - first index goes from 0 to N+1, 0=N, 1=N+1, the second index goes from
// 0 to m*o
// - array[i][j]= -{(\pi*(N-k)/N)/(m*o)}*j+\pi*(N-k)/N,   1<=i<=k
// - array[i][j]= -{(\pi*(-k)/N)/(m*o)}*j+\pi*(-k)/N,   1<=i<=k

// float J not needed for linear initial configuration but we have just included
// for fun.



#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "set_init.h"


void set_init(mpfr_t **array, int N,int k,int m,int o,mpfr_t J,int prec)
{
// variables

	mpfr_t pi;
	mpfr_t num1;
	mpfr_t num2;
	mpfr_t den;
	mpfr_t num1dden;
	mpfr_t num2dden;
	
//	mpfr_t **array = new mpfr_t * [N+2];
//	array[0] = new mpfr_t [(N+2)*(m*o+1)];
	
//	for (int i=1; i<N+2; i++){
//		array[i] = &array[0][i*(m*o+1)];
//	}
	
//	array[N] = &array[0][0];
//	array[N+1] = &array[0][m*o+1];
	
	mpfr_t *init0 = new mpfr_t [N+2];

	mpfr_t *initJ = new mpfr_t [N+2];
	
// variables
	
// initialization

//	for (int p=0; p<N+2; p++){
//		for (int l=0; l<m*o+1; l++){
//			mpfr_init2(array[p][l],prec);
//		}
//	}
	
	for (int p=0; p<N+2; p++){
		mpfr_init2(init0[p],prec);
		mpfr_init2(initJ[p],prec);
	}
	
	mpfr_init2(pi,prec);
	mpfr_init2(num1,prec);
	mpfr_init2(num2,prec);
	mpfr_init2(den,prec);
	mpfr_init2(num1dden,prec);
	mpfr_init2(num2dden,prec);

// initialization
	
// setting the initial value
		
	mpfr_const_pi(pi,MPFR_RNDD);
	
	for (int p=0; p<N+2; p++){
		mpfr_set_d(initJ[p],0.0,MPFR_RNDD);
	}	

	mpfr_set_si(num1,N-k,MPFR_RNDD);
	mpfr_set_si(num2,k,MPFR_RNDD);
	mpfr_neg(num2,num2,MPFR_RNDD);
	mpfr_set_si(den,N,MPFR_RNDD);
	mpfr_div(num1dden,num1,den,MPFR_RNDD);
	mpfr_div(num2dden,num2,den,MPFR_RNDD);

	
	for (int p=1; p<k+1; p++){
		mpfr_mul(init0[p],pi,num1dden,MPFR_RNDD);
	}
	for (int p=k+1; p<N+1; p++){
		mpfr_mul(init0[p],pi,num2dden,MPFR_RNDD);
	}
	
	for (int p=1; p<N+1; p++){
		mpfr_set(array[p][0],init0[p],MPFR_RNDD);
		mpfr_set(array[p][m*o],initJ[p],MPFR_RNDD);
		}

	for (int p=1; p<N+1; p++){
		for (int j=1; j<m*o; j++){
			mpfr_div_si(array[p][j],init0[p],m*o,MPFR_RNDD);
			mpfr_neg(array[p][j],array[p][j],MPFR_RNDD);
			mpfr_mul_si(array[p][j],array[p][j],j,MPFR_RNDD);
			mpfr_add(array[p][j],array[p][j],init0[p],MPFR_RNDD);
	}
	}
	
// setting the initial value

}
