// Calculates the contribution of the ij-th component of the array to the action
// of the string tension in making a first order approximation.

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "maction1ij.h"



void action1ij(mpfr_t** a,mpfr_t * ret, int N, int m, int o, int i, int j, mpfr_t J, mpfr_t zero, int prec)
{


// variables
	
//	mpfr_t* ret=new mpfr_t [1];
	mpfr_t S1;
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
	mpfr_init(S2);
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(num);
	mpfr_init(num1);
	mpfr_init(num2);
	mpfr_init(den);
	mpfr_init(ratio);
//	mpfr_init(ret[0]);
	mpfr_init(denabs);

// initialization

// calculating S1

	mpfr_set(S1,a[i][j],MPFR_RNDD);
	mpfr_sub(S1,S1,a[i][j+1],MPFR_RNDD);
	mpfr_sub(S1,S1,a[i][j-1],MPFR_RNDD);
	mpfr_mul(S1,S1,a[i][j],MPFR_RNDD);
	mpfr_mul_si(S1,S1,2*m,MPFR_RNDD);
	mpfr_div(S1,S1,J,MPFR_RNDD);
	
// calculating S1
	
// calculating S2

	mpfr_set_d(S2,4.0,MPFR_RNDD);

	for (int l=0;l<2;l++){
		for (int k=0;k<2;k++){
	
	mpfr_sub(arg1,a[i-l][j+k],a[i+1-l][j+k],MPFR_RNDD);
	mpfr_sub(arg2,a[i-l][j-1+k],a[i+1-l][j-1+k],MPFR_RNDD);

    mpfr_sub(den,arg1,arg2,MPFR_RNDD);

    mpfr_abs(denabs,den,MPFR_RNDD);
            
	if (mpfr_cmp(denabs,zero)<=0){
		mpfr_cos(num,arg2,MPFR_RNDD);
		mpfr_sub(S2,S2,num,MPFR_RNDD);
	}
	else {
		mpfr_sin(num1,arg1,MPFR_RNDD);
		mpfr_sin(num2,arg2,MPFR_RNDD);
		mpfr_sub(num,num1,num2,MPFR_RNDD);
		mpfr_div(ratio,num,den,MPFR_RNDD);
		mpfr_sub(S2,S2,ratio,MPFR_RNDD);
	}
		}
	}
	
    mpfr_mul(S2,S2,J,MPFR_RNDD);
	mpfr_div_ui(S2,S2,m,MPFR_RNDD);


// calculating S2
		
	mpfr_add(*ret,S1,S2,MPFR_RNDD);

//	mpfr_set(ret[0],S2,MPFR_RNDD);
    
// clearing variables
	
	mpfr_clear(S1);
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
	
}

