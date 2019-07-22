// Calculating the action to second order for the ij-th term. This version calculates the
// S2 term correctly as checked by the python code. S1 also matches.

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>
#include "fresnel.h"



void action2ij(mpfr_t** a, mpfr_t * ret, int N, int m, int o, int i, int j, mpfr_t J, mpfr_t** A, mpfr_t** M, mpfr_t zero, int prec)
{

// variables
    
    mpfr_t P1;
    mpfr_t P2;

	mpfr_t H1;
	mpfr_t H2;

    mpfr_t I;
    mpfr_t Dx;

	mpfr_t S1;
    mpfr_t S2;

    mpfr_t* w=new mpfr_t [3];

    mpfr_t comw1;
    mpfr_t comw2;

    mpfr_t sw2;
    
    mpfr_t arg1;
	mpfr_t arg2;
	mpfr_t arg3;
    
    int Numf=70;
	int r;
    
    char* ch=new char [1]; // debugging purposes
    
    mpfr_exp_t exp; // debugging purposes

// variables

// initialization
    
	mpfr_init(P1);
	mpfr_init(P2);
    
	mpfr_init(H1);
	mpfr_init(H2);

	mpfr_init(I);
	mpfr_init(Dx);
    
    mpfr_init(S1);
	mpfr_init(S2);

    for (int y=0;y<3;y++)
    {
        mpfr_init(w[y]);
    }
	
    mpfr_init(comw1);
    mpfr_init(comw2);
    
    mpfr_init(sw2);
    
	mpfr_init(arg1);
	mpfr_init(arg2);
	mpfr_init(arg3);
    
// initialization
    
    mpfr_set(Dx,J,MPFR_RNDD);
    mpfr_div_ui(Dx,Dx,m,MPFR_RNDD);
    r=j-(j/o)*o;

// calculating S1
    
    mpfr_set_flt(P1,0.0,MPFR_RNDD);
    mpfr_set_flt(S1,0.0,MPFR_RNDD);
    
    if (r==0)
    {
        for (int k=0; k<2; k++) {
            mpfr_set_flt(P1,0.0,MPFR_RNDD);
            for (int p=0; p<3; p++) {
                for (int q=0; q<3; q++) {
                    mpfr_mul(H1,a[i][j-k*o+p],a[i][j-k*o+q],MPFR_RNDD);
                    mpfr_mul(H1,H1,A[p][q],MPFR_RNDD);
                    mpfr_add(P1,P1,H1,MPFR_RNDD);
                }
            }
        mpfr_add(S1,S1,P1,MPFR_RNDD);
        }
        mpfr_div(S1,S1,Dx,MPFR_RNDD);
    }
    else if (r==1)
    {
        for (int p=0; p<3; p++) {
            for (int q=0; q<3; q++) {
                mpfr_mul(H1,a[i][j-r+p],a[i][j-r+q],MPFR_RNDD);
                mpfr_mul(H1,H1,A[p][q],MPFR_RNDD);
                mpfr_add(S1,S1,H1,MPFR_RNDD);
            }
        }
        mpfr_div(S1,S1,Dx,MPFR_RNDD);
    }
    
// calculating S1
    
	
// calculating S2
    
    mpfr_set_d(S2,0.0,MPFR_RNDD);

    if (r==0)
    {
        for (int q=0;q<2;q++)
        {
            for (int p=0;p<2;p++)
            {
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
//                        mpfr_get_str(ch,&exp,10,0,w[l],MPFR_RNDD);     //debugging purposes
                    }
                }
                mpfr_sub(H2,a[i-q][j-p*o],a[i-q+1][j-p*o],MPFR_RNDD);
                mpfr_mul(w[0],H2,M[0][0],MPFR_RNDD);
//                mpfr_get_str(ch,&exp,10,0,w[0],MPFR_RNDD);     //debugging purposes

                mpfr_abs(comw1,w[1],MPFR_RNDD);

//                mpfr_get_str(ch,&exp,10,0,comw1,MPFR_RNDD);     //debugging purposes
                
                mpfr_abs(comw2,w[2],MPFR_RNDD);

//                mpfr_get_str(ch,&exp,10,0,comw2,MPFR_RNDD);     //debugging purposes
                
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
    }
    else if (r==1)
    {
        for (int q=0;q<2;q++)
        {
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
//                        mpfr_get_str(ch,&exp,10,0,w[l],MPFR_RNDD);     //debugging purposes
                    }
                }
                mpfr_sub(H2,a[i-q][j-r],a[i-q+1][j-r],MPFR_RNDD);
                mpfr_mul(w[0],H2,M[0][0],MPFR_RNDD);
//                mpfr_get_str(ch,&exp,10,0,w[0],MPFR_RNDD);     //debugging purposes
            
                mpfr_abs(comw1,w[1],MPFR_RNDD);
            
//                mpfr_get_str(ch,&exp,10,0,comw1,MPFR_RNDD);     //debugging purposes

                mpfr_abs(comw2,w[2],MPFR_RNDD);

//                mpfr_get_str(ch,&exp,10,0,comw2,MPFR_RNDD);     //debugging purposes

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
		
    mpfr_add(*ret,S1,S2,MPFR_RNDD);

//    mpfr_set(*ret,S1,MPFR_RNDD);

// clearing variables
	
	mpfr_clear(P1);
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
        
    delete [] ch;
    
// clearing variables

}

