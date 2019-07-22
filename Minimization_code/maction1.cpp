// Calculates the string tension action to 1st order for an array with boundary
// conditions a[i][0]=\pi*(\mu_k)_i, a[i][m*o]=0. getting the array, number of gauge
// group, number of partition, order (its 1 but just for fun!), boundary integer J,
// precision of calculations and spits out the value of the action.

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "maction1.h"



mpfr_t* action1(mpfr_t** a, int N, int m, int o, mpfr_t J,mpfr_t zero, int prec)
{


// variables
	
	mpfr_t* ret=new mpfr_t [1];
	mpfr_t H1;
	mpfr_t S1;
	mpfr_t H2;
	mpfr_t S2;
	mpfr_t arg1;
	mpfr_t arg2;
	mpfr_t num;
	mpfr_t num1;
	mpfr_t num2;
	mpfr_t den;
	mpfr_t ratio;
	mpfr_t denabs;
	
// variables

// initialization
	
	mpfr_init(S1);
	mpfr_init(H1);
	mpfr_init(S2);
	mpfr_init(H2);
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(num);
	mpfr_init(num1);
	mpfr_init(num2);
	mpfr_init(den);
	mpfr_init(ratio);
	mpfr_init(ret[0]);
	mpfr_init(denabs);

// initialization

// calculating S1
	
	mpfr_set_d(S1,0.0,MPFR_RNDD);

	for (int i=1;i<N+1;i++){
		for (int j=1;j<m*o+1;j++){
			mpfr_sub(H1,a[i][j],a[i][j-1],MPFR_RNDD);
			mpfr_mul(H1,H1,H1,MPFR_RNDD);
			mpfr_mul_si(H1,H1,m,MPFR_RNDD);
			mpfr_div(H1,H1,J,MPFR_RNDD);
			mpfr_add(S1,S1,H1,MPFR_RNDD);
		}
	}
// calculating S1
	
// calculating S2

	mpfr_set_d(S2,0.0,MPFR_RNDD);

	for (int i=1;i<N+1;i++){
		for (int j=1;j<m*o+1;j++){
			mpfr_sub(arg1,a[i][j],a[i+1][j],MPFR_RNDD);
			mpfr_sub(arg2,a[i][j-1],a[i+1][j-1],MPFR_RNDD);
    mpfr_sub(den,arg1,arg2,MPFR_RNDD);
            
    mpfr_abs(denabs,den,MPFR_RNDD);

	if (mpfr_cmp(denabs,zero)<=0){
        mpfr_cos(num,arg2,MPFR_RNDD);
		mpfr_d_sub(H2,1.0,num,MPFR_RNDD);
	}
	else {
		mpfr_sin(num1,arg1,MPFR_RNDD);
		mpfr_sin(num2,arg2,MPFR_RNDD);
		mpfr_sub(num,num1,num2,MPFR_RNDD);
		mpfr_sub(den,arg1,arg2,MPFR_RNDD);
		mpfr_div(ratio,num,den,MPFR_RNDD);
		mpfr_d_sub(H2,1.0,ratio,MPFR_RNDD);
	}
		mpfr_mul(H2,H2,J,MPFR_RNDD);
		mpfr_div_ui(H2,H2,m,MPFR_RNDD);

		mpfr_add(S2,S2,H2,MPFR_RNDD);

		}
	}


// calculating S2
			
	mpfr_add(ret[0],S1,S2,MPFR_RNDD);

// clearing variables	
	
	mpfr_clear(H1);
	mpfr_clear(S1);
	mpfr_clear(H2);
	mpfr_clear(S2);
	mpfr_clear(arg1);
	mpfr_clear(arg2);
	mpfr_clear(num);
	mpfr_clear(num1);
	mpfr_clear(num2);
	mpfr_clear(den);
	mpfr_clear(ratio);
	mpfr_clear(denabs);
	
// clearing variables	
	
	return ret;
}

