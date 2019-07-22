// calcualting the fresnelcos integral

#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "fresnel.h"



void cos2int(mpfr_t* I, mpfr_t c2, mpfr_t c1, mpfr_t c0, mpfr_t zero, int prec)
{

// variables
    
    mpfr_t P;
    
	mpfr_t H;

    mpfr_t absa2;

    mpfr_t sa2;
    
    mpfr_t arg1;
	mpfr_t arg2;
	mpfr_t arg3;

    mpfr_t a0;
	mpfr_t a1;
	mpfr_t a2;
    
    
    
    int Numf=100;
    
    int z;
    int z2;
    
// variables

// initialization
    
	mpfr_init(P);
    
	mpfr_init(H);

    mpfr_init(sa2);
    
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(arg3);

	mpfr_init(a0);
	mpfr_init(a1);
	mpfr_init(a2);
    
    mpfr_init(absa2);
    
// initialization

    
    mpfr_set(a0,c0,MPFR_RNDD);
	mpfr_set(a1,c1,MPFR_RNDD);
	mpfr_set(a2,c2,MPFR_RNDD);


// calculating the fresnelcos integral
    
    mpfr_abs(absa2, a2, MPFR_RNDD);
    
    z=mpfr_cmp(absa2,zero);
    
    z2=mpfr_cmp(a2, zero);
    
    
    if (z<=0) {
        std::cout << "This module not intended for a2=0" << std::endl;
    }
    
    else if (z>0 and z2<0)
    {
    mpfr_neg(a2,a2,MPFR_RNDD);
    mpfr_neg(a1,a1,MPFR_RNDD);
    mpfr_neg(a0,a0,MPFR_RNDD);
        
    mpfr_sqrt(sa2,a2,MPFR_RNDD);

    //arg1
    mpfr_set(arg1,a0,MPFR_RNDD);
    mpfr_mul(H,a1,a1,MPFR_RNDD);
    mpfr_div(H,H,a2,MPFR_RNDD);
    mpfr_div_ui(H,H,4,MPFR_RNDD);
    mpfr_sub(arg1,arg1,H,MPFR_RNDD);
    //arg1
                    
    //arg2
    mpfr_set(arg2,sa2,MPFR_RNDD);
    mpfr_set(H,a1,MPFR_RNDD);
    mpfr_div(H,H,sa2,MPFR_RNDD);
    mpfr_div_ui(H,H,2,MPFR_RNDD);
    mpfr_add(arg2,arg2,H,MPFR_RNDD);
    //arg2
                    
    //arg3
    mpfr_set(arg3,a1,MPFR_RNDD);
    mpfr_div(arg3,arg3,sa2,MPFR_RNDD);
    mpfr_div_ui(arg3,arg3,2,MPFR_RNDD);
    //arg3
                    
    mpfr_cos(H, arg1, MPFR_RNDD);
    mpfr_sub(*I,C(arg2, Numf, prec)[0],C(arg3, Numf, prec)[0],MPFR_RNDD);
    mpfr_mul(*I,*I,H,MPFR_RNDD);
                    
    mpfr_sin(H, arg1, MPFR_RNDD);
    mpfr_sub(P,S(arg2, Numf, prec)[0],S(arg3, Numf, prec)[0],MPFR_RNDD);
    mpfr_mul(P,P,H,MPFR_RNDD);
                    
    mpfr_sub(*I,*I,P,MPFR_RNDD);
                    
    mpfr_div(*I,*I,sa2,MPFR_RNDD);

    }
    else if (z>0 and z2>0)
    {
        mpfr_sqrt(sa2,a2,MPFR_RNDD);
        
        //arg1
        mpfr_set(arg1,a0,MPFR_RNDD);
        mpfr_mul(H,a1,a1,MPFR_RNDD);
        mpfr_div(H,H,a2,MPFR_RNDD);
        mpfr_div_ui(H,H,4,MPFR_RNDD);
        mpfr_sub(arg1,arg1,H,MPFR_RNDD);
        //arg1
        
        //arg2
        mpfr_set(arg2,sa2,MPFR_RNDD);
        mpfr_set(H,a1,MPFR_RNDD);
        mpfr_div(H,H,sa2,MPFR_RNDD);
        mpfr_div_ui(H,H,2,MPFR_RNDD);
        mpfr_add(arg2,arg2,H,MPFR_RNDD);
        //arg2
        
        //arg3
        mpfr_set(arg3,a1,MPFR_RNDD);
        mpfr_div(arg3,arg3,sa2,MPFR_RNDD);
        mpfr_div_ui(arg3,arg3,2,MPFR_RNDD);
        //arg3
        
        mpfr_cos(H, arg1, MPFR_RNDD);
        mpfr_sub(*I,C(arg2, Numf, prec)[0],C(arg3, Numf, prec)[0],MPFR_RNDD);
        mpfr_mul(*I,*I,H,MPFR_RNDD);
        
        mpfr_sin(H, arg1, MPFR_RNDD);
        mpfr_sub(P,S(arg2, Numf, prec)[0],S(arg3, Numf, prec)[0],MPFR_RNDD);
        mpfr_mul(P,P,H,MPFR_RNDD);
        
        mpfr_sub(*I,*I,P,MPFR_RNDD);
        
        mpfr_div(*I,*I,sa2,MPFR_RNDD);
    }

// calculating the fresnelcos integral


// clearing variables
	
	mpfr_clear(P);

	mpfr_clear(H);

	mpfr_clear(sa2);

	mpfr_clear(arg1);
	mpfr_clear(arg2);
	mpfr_clear(arg3);

	mpfr_clear(a0);
	mpfr_clear(a1);
	mpfr_clear(a2);
    
// clearing variables

}

