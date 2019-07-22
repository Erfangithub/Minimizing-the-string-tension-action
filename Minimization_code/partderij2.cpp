// calculates the first partial derivative of the discretized action to second order  with respect to the ij-th variable

#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "fresnel.h"
#include "mfintz2sin.h"
#include "mfintzsin.h"
#include "fsinint.h"



void parder2ij(mpfr_t** a, mpfr_t * par, int N, int m, int o, int i, int j, mpfr_t J,mpfr_t** M, mpfr_t zero, int prec)
{

// variables
    
    mpfr_t P1;
    mpfr_t P2;

	mpfr_t H1;
	mpfr_t H2;
    
    mpfr_t csum;
    mpfr_t ssum;
    
    mpfr_t Dx;

	mpfr_t S1;
    mpfr_t S2;

    mpfr_t* w=new mpfr_t [3];

    mpfr_t comw1;
    mpfr_t comw2;

	int r;

// variables

// initialization
    
	mpfr_init(P1);
	mpfr_init(P2);
    
	mpfr_init(H1);
	mpfr_init(H2);

	mpfr_init(Dx);
    
    mpfr_init(S1);
	mpfr_init(S2);
    
	mpfr_init(csum);
	mpfr_init(ssum);
    

    for (int y=0;y<3;y++)
    {
        mpfr_init(w[y]);
    }
	
    mpfr_init(comw1);
    mpfr_init(comw2);
    
// initialization
    
    mpfr_set(Dx,J,MPFR_RNDD);
    mpfr_div_ui(Dx,Dx,m,MPFR_RNDD);
    r=j-(j/o)*o;

// calculating dijS1
        
    if (r==0)
    {
        mpfr_mul_ui(H1,a[i][j],28,MPFR_RNDD);

        mpfr_add(P1,a[i][j+1],a[i][j-1],MPFR_RNDD);
        
        mpfr_mul_ui(P1,P1,16,MPFR_RNDD);
        
        mpfr_sub(H1,H1,P1,MPFR_RNDD);

        mpfr_add(P1,a[i][j+2],a[i][j-2],MPFR_RNDD);

        mpfr_mul_ui(P1,P1,2,MPFR_RNDD);
        
        mpfr_add(H1,H1,P1,MPFR_RNDD);
        
        mpfr_div_ui(H1,H1,3,MPFR_RNDD);

        mpfr_div(S1,H1,Dx,MPFR_RNDD);
    }
    else if (r==1)
    {
        mpfr_mul_ui(H1,a[i][j],32,MPFR_RNDD);
        
        mpfr_add(P1,a[i][j+1],a[i][j-1],MPFR_RNDD);
        
        mpfr_mul_ui(P1,P1,16,MPFR_RNDD);
        
        mpfr_sub(H1,H1,P1,MPFR_RNDD);

        mpfr_div_ui(H1,H1,3,MPFR_RNDD);
        
        mpfr_div(S1,H1,Dx,MPFR_RNDD);
    }
    
// calculating dijS1
    
	
// calculating dijS2
    
    mpfr_set_d(S2,0.0,MPFR_RNDD);

    if (r==0)
    {
        for (int q=0;q<2;q++) // i=q
        {
            for (int p=0;p<2;p++) // j=p
            {
                // calculating the w coefficients

                for (int y=0;y<3;y++)
                {
                    mpfr_set_d(w[y],0.0,MPFR_RNDD);
                }

                for (int l=1;l<3;l++)
                {
                    for (int k=0;k<3;k++)
                    {
                        mpfr_sub(H2,a[i-q][j-p*o+k],a[i-q+1][j-p*o+k],MPFR_RNDD);
                        mpfr_mul(H2,H2,M[l][k],MPFR_RNDD);
                        mpfr_add(w[l],w[l],H2,MPFR_RNDD);
                    }
                }
                mpfr_sub(H2,a[i-q][j-p*o],a[i-q+1][j-p*o],MPFR_RNDD);
                mpfr_mul(w[0],H2,M[0][0],MPFR_RNDD);

                // calculating the w coefficients
                
                mpfr_abs(comw1,w[1],MPFR_RNDD);
                
                mpfr_abs(comw2,w[2],MPFR_RNDD);

                
                int z=mpfr_cmp(comw2,zero);

                int z2=mpfr_cmp(w[2],zero);
                
                if (z<=0)
                {
                    if (mpfr_cmp(comw1,zero)<=0)
                    {
                        mpfr_sin(P2,w[0],MPFR_RNDD);
                        mpfr_div_ui(P2,P2,6,MPFR_RNDD);
                        mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                        if (q==0) {
                            mpfr_add(S2,S2,P2,MPFR_RNDD);
                        }
                        else
                        {
                            mpfr_sub(S2,S2,P2,MPFR_RNDD);
                        }
                    }
                    else
                    {
                        // calculating int_0^1dzz^2sin(w1z+w0)
                        
                        mpfr_add(H2,w[1],w[0],MPFR_RNDD);
 
                        mpfr_cos(csum,H2,MPFR_RNDD);

                        mpfr_sin(ssum,H2,MPFR_RNDD);
                        
                        mpfr_cos(H2,w[0],MPFR_RNDD);
                        mpfr_sub(H2,csum,H2,MPFR_RNDD);
                        mpfr_div(H2,H2,w[1],MPFR_RNDD);
                        
                        mpfr_add(H2,H2,ssum,MPFR_RNDD);
                        mpfr_div(H2,H2,w[1],MPFR_RNDD);
                        mpfr_mul_ui(H2,H2,2,MPFR_RNDD);
                        mpfr_sub(H2,H2,csum,MPFR_RNDD);
                        
                        mpfr_div(H2,H2,w[1],MPFR_RNDD);

                        // calculating int_0^1dzz^2sin(w1z+w0)
                        
                        mpfr_mul_ui(P2,H2,2,MPFR_RNDD);
                        
                        // calculating int_0^1dzzsin(w1z+w0)
                        
                        mpfr_sin(H2,w[0],MPFR_RNDD);
                        mpfr_sub(H2,ssum,H2,MPFR_RNDD);
                        mpfr_div(H2,H2,w[1],MPFR_RNDD);
                        
                        mpfr_sub(H2,H2,csum,MPFR_RNDD);
                        mpfr_div(H2,H2,w[1],MPFR_RNDD);
                        
                        // calculating int_0^1dzzsin(w1z+w0)
                        
                        if (p==0) {
                            mpfr_mul_ui(H2,H2,3,MPFR_RNDD);
                            mpfr_sub(P2,P2,H2,MPFR_RNDD);
                            mpfr_cos(H2,w[0],MPFR_RNDD);
                            mpfr_sub(H2,H2,csum,MPFR_RNDD);
                            mpfr_div(H2,H2,w[1],MPFR_RNDD);
                            mpfr_add(P2,P2,H2,MPFR_RNDD);
                        }
                        else
                        {
                            mpfr_sub(P2,P2,H2,MPFR_RNDD);
                        }
                        mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                        
                        if (q==0) {
                            mpfr_add(S2,S2,P2,MPFR_RNDD);
                        }
                        else
                        {
                            mpfr_sub(S2,S2,P2,MPFR_RNDD);
                        }
                    }
                }
                else if (z>0)
                {
                    if (p==0) {
                        
                    fintz2sin(&H2, w[2], w[1], w[0], zero, prec);
                    
                    mpfr_mul_ui(P2,H2,2,MPFR_RNDD);

                    fintzsin(&H2, w[2], w[1], w[0], zero, prec);
                    
                    mpfr_mul_ui(H2,H2,3,MPFR_RNDD);

                    mpfr_sub(P2,P2,H2,MPFR_RNDD);

                    sin2int(&H2, w[2], w[1], w[0], zero, prec);

                    mpfr_add(P2,P2,H2,MPFR_RNDD);
                    
                    mpfr_mul(P2,P2,Dx,MPFR_RNDD);

                    }
                    else
                    {
                        fintz2sin(&H2, w[2], w[1], w[0], zero, prec);
                        
                        mpfr_mul_ui(P2,H2,2,MPFR_RNDD);
                        
                        fintzsin(&H2, w[2], w[1], w[0], zero, prec);

                        mpfr_sub(P2,P2,H2,MPFR_RNDD);
                        
                        mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    }
                    
                    if (q==0) {
                        mpfr_add(S2,S2,P2,MPFR_RNDD);
                    }
                    else
                    {
                        mpfr_sub(S2,S2,P2,MPFR_RNDD);
                    }
                }

            }
        }
    }
    else if (r==1)
    {
        for (int q=0;q<2;q++)
        {
            //calculating the w coefficients
            for (int y=0;y<3;y++)
            {
                mpfr_set_d(w[y],0.0,MPFR_RNDD);
            }
                
            for (int l=1;l<3;l++)
            {
                for (int k=0;k<3;k++)
                {
                    mpfr_sub(H2,a[i-q][j-r+k],a[i-q+1][j-r+k],MPFR_RNDD);
                    mpfr_mul(H2,H2,M[l][k],MPFR_RNDD);
                    mpfr_add(w[l],w[l],H2,MPFR_RNDD);
                }
            }
            mpfr_sub(H2,a[i-q][j-r],a[i-q+1][j-r],MPFR_RNDD);
            mpfr_mul(w[0],H2,M[0][0],MPFR_RNDD);
            
            //calculating the w coefficients

            mpfr_abs(comw1,w[1],MPFR_RNDD);
            
            mpfr_abs(comw2,w[2],MPFR_RNDD);

            int z=mpfr_cmp(comw2,zero);
            int z2=mpfr_cmp(w[2],zero);

            if (z<=0)
            {
                if (mpfr_cmp(comw1,zero)<=0)
                {
                    mpfr_sin(P2,w[0],MPFR_RNDD);
                    mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    mpfr_mul_ui(P2,P2,2,MPFR_RNDD);
                    mpfr_div_ui(P2,P2,3,MPFR_RNDD);
                    if (q==0) {
                        mpfr_add(S2,S2,P2,MPFR_RNDD);
                    }
                    else
                    {
                        mpfr_sub(S2,S2,P2,MPFR_RNDD);
                    }
                }
                else
                {
                    // calculating int_0^1dzz^2sin(w1z+w0)
                    
                    mpfr_add(H2,w[1],w[0],MPFR_RNDD);
                    
                    mpfr_cos(csum,H2,MPFR_RNDD);
                    
                    mpfr_sin(ssum,H2,MPFR_RNDD);
                    
                    mpfr_cos(H2,w[0],MPFR_RNDD);
                    mpfr_sub(H2,csum,H2,MPFR_RNDD);
                    mpfr_div(H2,H2,w[1],MPFR_RNDD);
                    
                    mpfr_add(H2,H2,ssum,MPFR_RNDD);
                    mpfr_div(H2,H2,w[1],MPFR_RNDD);
                    mpfr_mul_ui(H2,H2,2,MPFR_RNDD);
                    mpfr_sub(H2,H2,csum,MPFR_RNDD);
                    
                    mpfr_div(P2,H2,w[1],MPFR_RNDD);
                    
                    // calculating int_0^1dzz^2sin(w1z+w0)
                                        
                    // calculating int_0^1dzzsin(w1z+w0)
                    
                    mpfr_sin(H2,w[0],MPFR_RNDD);
                    mpfr_sub(H2,ssum,H2,MPFR_RNDD);
                    mpfr_div(H2,H2,w[1],MPFR_RNDD);
                    
                    mpfr_sub(H2,H2,csum,MPFR_RNDD);
                    mpfr_div(H2,H2,w[1],MPFR_RNDD);
                    
                    // calculating int_0^1dzzsin(w1z+w0)
                    
                    mpfr_sub(P2,H2,P2,MPFR_RNDD);
                    mpfr_mul_ui(P2,P2,4,MPFR_RNDD);
                    mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                    
                    if (q==0) {
                        mpfr_add(S2,S2,P2,MPFR_RNDD);
                    }
                    else
                    {
                        mpfr_sub(S2,S2,P2,MPFR_RNDD);
                    }

                }
            }
            else if (z>0)
            {
                fintz2sin(&H2, w[2], w[1], w[0], zero, prec);
                fintzsin(&P2, w[2], w[1], w[0], zero, prec);
                mpfr_sub(P2,P2,H2,MPFR_RNDD);
                mpfr_mul_ui(P2,P2,4,MPFR_RNDD);
                mpfr_mul(P2,P2,Dx,MPFR_RNDD);
                if (q==0) {
                    mpfr_add(S2,S2,P2,MPFR_RNDD);
                }
                else
                {
                    mpfr_sub(S2,S2,P2,MPFR_RNDD);
                }
            }
            
        }
    }
    
// calculating dijS2
		
//    mpfr_add(*par,S1,S2,MPFR_RNDD);

    mpfr_set(*par,S2,MPFR_RNDD);
    

//    mpfr_set(*par,w[1],MPFR_RNDD);

// clearing variables
	
	mpfr_clear(P1);
	mpfr_clear(P2);

	mpfr_clear(H1);
	mpfr_clear(H2);

	mpfr_clear(Dx);

    mpfr_clear(S1);
	mpfr_clear(S2);

    mpfr_clear(csum);
	mpfr_clear(ssum);
    
    mpfr_clear(w[0]);
    mpfr_clear(w[1]);
    mpfr_clear(w[2]);
    
    mpfr_clear(comw1);
	mpfr_clear(comw2);
            
// clearing variables

}

