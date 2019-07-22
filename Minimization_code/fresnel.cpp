// Arbitrary precision Fresnel integrals:
//  C(x)=\int_0^xcos((\pi/2)t^2)dt
//  S(x)=\int_0^xsin((\pi/2)t^2)dt
//  gets x, number of terms in the series expansion of the above integrals,
//  precision of calculations and gives out C(x) and S(x).


#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>


// cos fresnel function

mpfr_t* C(mpfr_t x,int n, int prec)
{
// variables
	
	mpfr_t fac;
	mpfr_t value;
	mpfr_t var;
	mpfr_t *array=new mpfr_t [1];
	mpfr_t pi;
	mpfr_t coeff;
	mpfr_t xc;
	
// variables

// initialization

	mpfr_init(fac);
	mpfr_init(value);
	mpfr_init(var);
	mpfr_init(array[0]);
	mpfr_init(pi);
	mpfr_init(coeff);
	mpfr_init(xc);
	
// initialization

// summing the fresnel series
	
	mpfr_const_pi(pi,MPFR_RNDD);
	
	mpfr_div_d(coeff,pi,2.0,MPFR_RNDD);

	mpfr_sqrt(coeff,coeff,MPFR_RNDD);
	
	mpfr_mul(xc,x,coeff,MPFR_RNDD);
	
	mpfr_set_d(value,0.0,MPFR_RNDD);
		

	for (unsigned long int r=0;r<n+1;r++)
{
		
//		mpfr_pow_ui(var,xc,4*r+1,MPFR_RNDD);     // uncomment for the \pi/2 version
        mpfr_pow_ui(var,x,4*r+1,MPFR_RNDD);
        mpfr_div_ui(var,var,4*r+1,MPFR_RNDD);
		mpfr_fac_ui (fac, 2*r,MPFR_RNDD);
		mpfr_div(var,var,fac,MPFR_RNDD);

		if (r-2*(r/2)>0.5)
		{
			mpfr_neg(var,var,MPFR_RNDD);
		}
		
		mpfr_add(value,value,var,MPFR_RNDD);
		
	}
//	mpfr_div(value,value,coeff,MPFR_RNDD);       // uncomment for the \pi/2 version

// summing the fresnel series

// clearing variables	
	
	mpfr_clear(fac);
	mpfr_clear(var);
	mpfr_set(array[0],value,MPFR_RNDD);
	mpfr_clear(value);
	mpfr_clear(pi);
	mpfr_clear(coeff);
	mpfr_clear(xc);

// clearing variables	
	
	return array;
}

// cos fresnel function




mpfr_t* S(mpfr_t x,int n, int prec)
{
// variables
	
	mpfr_t fac;
	mpfr_t value;
	mpfr_t var;
	mpfr_t *array1=new mpfr_t [1];
	mpfr_t pi;
	mpfr_t coeff;
	mpfr_t xc;
	
// variables
	
// initialization
	
	mpfr_init2(fac,prec);
	mpfr_init2(value,prec);
	mpfr_init2(var,prec);
	mpfr_init2(array1[0], prec);
	mpfr_init2(pi, prec);
	mpfr_init2(coeff, prec);
	mpfr_init2(xc, prec);
	
// initialization
	
// summing the fresnel series
	
	mpfr_const_pi(pi,MPFR_RNDD);
	
	mpfr_div_d(coeff,pi,2.0,MPFR_RNDD);
	
	mpfr_sqrt(coeff,coeff,MPFR_RNDD);
	
	mpfr_mul(xc,x,coeff,MPFR_RNDD);
	
	mpfr_set_d(value,0.0,MPFR_RNDD);
	
	
	for (unsigned long int r=0;r<n+1;r++)
	{
		
//		mpfr_pow_ui(var,xc,4*r+3,MPFR_RNDD);    // uncomment for the \pi/2 version
		mpfr_pow_ui(var,x,4*r+3,MPFR_RNDD);
		mpfr_div_ui(var,var,4*r+3,MPFR_RNDD);
		mpfr_fac_ui (fac, 2*r+1, MPFR_RNDD);
		mpfr_div(var,var,fac,MPFR_RNDD);
		
		if (r-2*(r/2)>0.5)
		{
			mpfr_neg(var,var,MPFR_RNDD);
		}
		
		mpfr_add(value,value,var,MPFR_RNDD);
		
	}

//	mpfr_div(value,value,coeff,MPFR_RNDD);    // uncomment for the \pi/2 version
	
// summing the fresnel series

// clearing variables
	
	mpfr_clear(fac);
	mpfr_clear(var);
	mpfr_set(array1[0],value,MPFR_RNDD);
	mpfr_clear(value);
	mpfr_clear(pi);
	mpfr_clear(coeff);
	mpfr_clear(xc);

// clearing variables	
	
	return array1;
}

// sin fresnel function


//int main()
//{
	
//	mpfr_t r;
	
//	mpfr_init2(r,200);

//	mpfr_set_d(r,2.0,MPFR_RNDD);
	
//	mpfr_out_str (stdout, 10, 0, C(r,100,600)[0], MPFR_RNDD);

//	putchar ('\n');
	
//	mpfr_out_str (stdout, 10, 0, S(r,100,600)[0], MPFR_RNDD);

//	putchar ('\n');
	
//	return 0;
//}
