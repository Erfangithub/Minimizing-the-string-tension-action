// calcualting the fresnelcos integral

#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "fresnel.h"



void fintzsin(mpfr_t* I, mpfr_t c2, mpfr_t c1, mpfr_t c0, mpfr_t zero, int prec)
{
    
    // variables
    
    mpfr_t P;
    
	mpfr_t H;
    
	mpfr_t K;

	mpfr_t a2;

	mpfr_t a1;
    
	mpfr_t a0;
    
    mpfr_t absa2;
    
    mpfr_t sa2;
    
    mpfr_t arg1;
	mpfr_t arg2;
	mpfr_t arg3;
    
    int Numf=100;
    
    int z;
    int z2;
    
    // variables
    
    // initialization
    
	mpfr_init(P);
    
	mpfr_init(H);
    
    mpfr_init(K);

    mpfr_init(a2);

    mpfr_init(a1);

    mpfr_init(a0);
    
    mpfr_init(sa2);
    
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(arg3);
    
    mpfr_init(absa2);
    
    // initialization
    
    mpfr_set(a2,c2,MPFR_RNDD);

    mpfr_set(a1,c1,MPFR_RNDD);

    mpfr_set(a0,c0,MPFR_RNDD);
    
    
    // calculating the fresnelz2sin integral
    
    mpfr_abs(absa2, a2, MPFR_RNDD);
    
    z=mpfr_cmp(absa2,zero);
    
    z2=mpfr_cmp(a2, zero);
    
    if (z<=0) {
        std::cout << "This module not intended for a2=0" << std::endl;
    }
    
    else
    {
        if (z>0 and z2<0){
            mpfr_neg(a2,a2,MPFR_RNDD);
            mpfr_neg(a1,a1,MPFR_RNDD);
            mpfr_neg(a0,a0,MPFR_RNDD);
        }
        
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
        
        // calculating I_{1}
        
        
        mpfr_sub(H,S(arg2, Numf, prec)[0],S(arg3, Numf, prec)[0],MPFR_RNDD);

        mpfr_mul(H,H,arg3,MPFR_RNDD);

        mpfr_pow_ui(P,arg3,2,MPFR_RNDD);
        
        mpfr_add(P,P,a2,MPFR_RNDD);

        mpfr_add(P,P,a1,MPFR_RNDD);
        mpfr_cos(P,P,MPFR_RNDD);
        mpfr_div_ui(P,P,2,MPFR_RNDD);

        mpfr_add(H,H,P,MPFR_RNDD);
        
        mpfr_pow_ui(P,arg3,2,MPFR_RNDD);
        mpfr_cos(P,P,MPFR_RNDD);

        mpfr_div_ui(P,P,2,MPFR_RNDD);

        mpfr_sub(H,H,P,MPFR_RNDD);

        mpfr_cos(P,arg1,MPFR_RNDD);
        
        mpfr_mul(H,H,P,MPFR_RNDD);

        mpfr_div(H,H,a2,MPFR_RNDD);
        
        mpfr_neg(*I,H,MPFR_RNDD);
        
        // calculating I_{1}
        
                
        // calculating I_{2}
        
        
        mpfr_sub(K,C(arg3, Numf, prec)[0],C(arg2, Numf, prec)[0],MPFR_RNDD);
        
        mpfr_mul(K,K,arg3,MPFR_RNDD);
        
        mpfr_pow_ui(P,arg3,2,MPFR_RNDD);
        
        mpfr_add(P,P,a2,MPFR_RNDD);
        
        mpfr_add(P,P,a1,MPFR_RNDD);
        mpfr_sin(P,P,MPFR_RNDD);
        mpfr_div_ui(P,P,2,MPFR_RNDD);
        
        mpfr_add(K,K,P,MPFR_RNDD);
        
        mpfr_pow_ui(P,arg3,2,MPFR_RNDD);
        mpfr_sin(P,P,MPFR_RNDD);
        
        mpfr_div_ui(P,P,2,MPFR_RNDD);
        
        mpfr_sub(K,K,P,MPFR_RNDD);
        
        mpfr_sin(P,arg1,MPFR_RNDD);
        
        mpfr_mul(K,K,P,MPFR_RNDD);
        
        mpfr_div(K,K,a2,MPFR_RNDD);
                
        // calculating I_{2}
        
        mpfr_add(*I,*I,K,MPFR_RNDD);

        if (z>0 and z2<0){
            mpfr_neg(*I,*I,MPFR_RNDD);
        }
    }
    
    // calculating the fresnelz2sin integral
    
    
    // clearing variables
	
	mpfr_clear(P);
    
	mpfr_clear(H);
    
	mpfr_clear(sa2);
    
	mpfr_clear(arg1);
	mpfr_clear(arg2);
	mpfr_clear(arg3);
    
    // clearing variables
    
}

