// Calculates the string tension action to 2nd order for an array with boundary
// conditions a[i][0]=\pi*(\mu_k)_i, a[i][m*o]=0. getting the array, number of gauge
// group, number of partition, order (its 2 but just for fun!), boundary integer J,
// precision of calculations and spits out the value of the action.

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "maction2.h"
#include "fresnel.h"



void action2(mpfr_t** a, mpfr_t* action , int N, int m, int o, mpfr_t J, mpfr_t** A, mpfr_t** M, mpfr_t zero, int prec)
{


// variables
	
	mpfr_t H1;
	mpfr_t S1;
	mpfr_t Dx;

	
    mpfr_t H2;
	mpfr_t S2;
    mpfr_t P2;
    mpfr_t I;
    mpfr_t arg1;
	mpfr_t arg2;
	mpfr_t arg3;
    
    mpfr_t* w=new mpfr_t [3];
    
    mpfr_t comw1;
    mpfr_t comw2;
    
    mpfr_t sw2;

    int Numf=100;
	
// variables

// initialization
	

	mpfr_init(S1);
	mpfr_init(H1);
	mpfr_init(Dx);
    
    
	mpfr_init(S2);
	mpfr_init(H2);
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(arg3);
	mpfr_init(P2);
	mpfr_init(I);
	mpfr_init(comw1);
	mpfr_init(comw2);
	mpfr_init(sw2);
    for (int y=0;y<3;y++)
    {
        mpfr_init(w[y]);
    }

// initialization
    
    mpfr_set(Dx,J,MPFR_RNDD);
    mpfr_div_ui(Dx,Dx,m,MPFR_RNDD);

// calculating S1
	
	mpfr_set_d(S1,0.0,MPFR_RNDD);

	for (int i=1;i<N+1;i++){
		for (int j=0;j<m;j++){
            for (int p=0; p<3; p++) {
                for (int q=0; q<3; q++) {
                    mpfr_mul(H1,a[i][j*o+p],a[i][j*o+q],MPFR_RNDD);
                    mpfr_mul(H1,H1,A[p][q],MPFR_RNDD);
                    mpfr_add(S1,S1,H1,MPFR_RNDD);
                }
            }
		}
	}
    mpfr_div(S1,S1,Dx,MPFR_RNDD);

// calculating S1
	
// calculating S2

    mpfr_set_flt(S2,0.0,MPFR_RNDD);
    
    for (int i=1;i<N+1;i++)
        {
            for (int j=0;j<m;j++)
            {
                for (int y=0;y<3;y++)
                {
                    mpfr_set_d(w[y],0.0,MPFR_RNDD);
                }
                
                for (int l=1;l<3;l++)
                {
                    for (int k=0;k<3;k++)
                    {
                        mpfr_sub(H2,a[i][j*o+k],a[i+1][j*o+k],MPFR_RNDD);
                        mpfr_mul(H2,H2,M[l][k],MPFR_RNDD);
                        mpfr_add(w[l],w[l],H2,MPFR_RNDD);
                    }
                }
                mpfr_sub(H2,a[i][j*o],a[i+1][j*o],MPFR_RNDD);
                mpfr_mul(w[0],H2,M[0][0],MPFR_RNDD);
                
                mpfr_abs(comw1,w[1],MPFR_RNDD);
                mpfr_abs(comw2,w[2],MPFR_RNDD);
                int z=mpfr_cmp(comw2,zero);
                int z2=mpfr_cmp(w[2],zero);
                
                if (z<=0)
                {
                    if (mpfr_cmp(comw1,zero)<=0)
                    {
                        mpfr_cos(P2,w[0],MPFR_RNDD);
                        mpfr_ui_sub(P2,1,P2,MPFR_RNDD);
                        mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    }
                    else
                    {
                        mpfr_add(H2,w[1],w[0],MPFR_RNDD);
                        mpfr_sin(P2,H2,MPFR_RNDD);
                        mpfr_set(H2,w[0],MPFR_RNDD);
                        mpfr_sin(H2,H2,MPFR_RNDD);
                        mpfr_sub(P2,P2,H2,MPFR_RNDD);
                        mpfr_div(P2,P2,w[1],MPFR_RNDD);
                        mpfr_ui_sub(P2,1,P2,MPFR_RNDD);
                        mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    }
                }
                else if (z>0 and z2<0)
                {
                    mpfr_neg(w[2],w[2],MPFR_RNDD);
                    mpfr_neg(w[1],w[1],MPFR_RNDD);
                    mpfr_neg(w[0],w[0],MPFR_RNDD);
                    
                    mpfr_sqrt(sw2,w[2],MPFR_RNDD);
                    //arg1
                    mpfr_set(arg1,w[0],MPFR_RNDD);
                    mpfr_mul(H2,w[1],w[1],MPFR_RNDD);
                    mpfr_div(H2,H2,w[2],MPFR_RNDD);
                    mpfr_div_ui(H2,H2,4,MPFR_RNDD);
                    mpfr_sub(arg1,arg1,H2,MPFR_RNDD);
                    //arg1
                    
                    //arg2
                    mpfr_set(arg2,sw2,MPFR_RNDD);
                    mpfr_set(H2,w[1],MPFR_RNDD);
                    mpfr_div(H2,H2,sw2,MPFR_RNDD);
                    mpfr_div_ui(H2,H2,2,MPFR_RNDD);
                    mpfr_add(arg2,arg2,H2,MPFR_RNDD);
                    //arg2
                    
                    //arg3
                    mpfr_set(arg3,w[1],MPFR_RNDD);
                    mpfr_div(arg3,arg3,sw2,MPFR_RNDD);
                    mpfr_div_ui(arg3,arg3,2,MPFR_RNDD);
                    //arg3
                    
                    mpfr_cos(H2, arg1, MPFR_RNDD);
                    mpfr_sub(I,C(arg2, Numf, prec)[0],C(arg3, Numf, prec)[0],MPFR_RNDD);
                    mpfr_mul(P2,I,H2,MPFR_RNDD);
                    
                    mpfr_sin(H2, arg1, MPFR_RNDD);
                    mpfr_sub(I,S(arg2, Numf, prec)[0],S(arg3, Numf, prec)[0],MPFR_RNDD);
                    mpfr_mul(I,I,H2,MPFR_RNDD);
                    
                    mpfr_sub(P2,P2,I,MPFR_RNDD);
                    
                    mpfr_div(P2,P2,sw2,MPFR_RNDD);
                    
                    mpfr_ui_sub(P2,1,P2,MPFR_RNDD);
                    
                    mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    
                }
                else if (z>0 and z2>0)
                {
                    mpfr_sqrt(sw2,w[2],MPFR_RNDD);
                    //arg1
                    mpfr_set(arg1,w[0],MPFR_RNDD);
                    mpfr_mul(H2,w[1],w[1],MPFR_RNDD);
                    mpfr_div(H2,H2,w[2],MPFR_RNDD);
                    mpfr_div_ui(H2,H2,4,MPFR_RNDD);
                    mpfr_sub(arg1,arg1,H2,MPFR_RNDD);
                    //arg1
                    
                    //arg2
                    mpfr_set(arg2,sw2,MPFR_RNDD);
                    mpfr_set(H2,w[1],MPFR_RNDD);
                    mpfr_div(H2,H2,sw2,MPFR_RNDD);
                    mpfr_div_ui(H2,H2,2,MPFR_RNDD);
                    mpfr_add(arg2,arg2,H2,MPFR_RNDD);
                    //arg2
                    
                    //arg3
                    mpfr_set(arg3,w[1],MPFR_RNDD);
                    mpfr_div(arg3,arg3,sw2,MPFR_RNDD);
                    mpfr_div_ui(arg3,arg3,2,MPFR_RNDD);
                    //arg3
                    
                    mpfr_cos(H2, arg1, MPFR_RNDD);
                    mpfr_sub(I,C(arg2, Numf, prec)[0],C(arg3, Numf, prec)[0],MPFR_RNDD);
                    mpfr_mul(P2,I,H2,MPFR_RNDD);
                    
                    mpfr_sin(H2, arg1, MPFR_RNDD);
                    mpfr_sub(I,S(arg2, Numf, prec)[0],S(arg3, Numf, prec)[0],MPFR_RNDD);
                    mpfr_mul(I,I,H2,MPFR_RNDD);
                    
                    mpfr_sub(P2,P2,I,MPFR_RNDD);
                    
                    mpfr_div(P2,P2,sw2,MPFR_RNDD);
                    
                    mpfr_ui_sub(P2,1,P2,MPFR_RNDD);
                    
                    mpfr_mul(P2,P2,Dx,MPFR_RNDD);

                }
                mpfr_add(S2,S2,P2,MPFR_RNDD);
            }
        }

// calculating S2
			
	mpfr_add(*action,S1,S2,MPFR_RNDD);

// clearing variables
	
	mpfr_clear(P2);
    
	mpfr_clear(H1);
	mpfr_clear(H2);
    
	mpfr_clear(I);
	mpfr_clear(Dx);
    
    mpfr_clear(S1);
	mpfr_clear(S2);
    
    mpfr_clear(w[0]);
    mpfr_clear(w[1]);
    mpfr_clear(w[2]);
    
    mpfr_clear(comw1);
	mpfr_clear(comw2);
    
	mpfr_clear(sw2);
    
	mpfr_clear(arg1);
	mpfr_clear(arg2);
	mpfr_clear(arg3);
    
// clearing variables
	
}

