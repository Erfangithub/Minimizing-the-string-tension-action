// second order version of main



#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include "set_init.h"
#include "maction1ij.h"
#include "maction1.h"
#include "maction2ij.h"
#include "maction2.h"
#include "partderij2.h"


int main(int argc, const char * argv[])
{
	
    // variables
    
    int prec=500;
    
    mpfr_set_default_prec (prec);
    
	int N=2;
	int k=1;
	int m=10;
	int o=2;
    
	mpfr_t J;
    
    mpfr_t par;
    
    mpfr_init(J);
    
    mpfr_set_ui(J,6,MPFR_RNDD);
    
    int numvar=N*(m*o-1);
    int maxj=m*o-1;
    
    mpfr_t action;
    
    mpfr_t c;
    
    mpfr_t zero;
    
    mpfr_t w;
    
    mpfr_t actija;
    
    mpfr_t actijt;
    
    mpfr_t **A = new mpfr_t * [o+1];
	A[0] = new mpfr_t [(o+1)*(o+1)];
    for (int i=1; i<o+1; i++){
		A[i] = &A[0][i*(o+1)];
	}
    
    mpfr_t **M = new mpfr_t * [o+1];
	M[0] = new mpfr_t [(o+1)*(o+1)];
    for (int i=1; i<o+1; i++){
		M[i] = &M[0][i*(o+1)];
    }
    
    // creating arrays
    
	mpfr_t **a = new mpfr_t * [N+2];
	a[0] = new mpfr_t [(N+2)*(m*o+1)];
	
	for (int i=1; i<N+2; i++){
		a[i] = &a[0][i*(m*o+1)];
	}
	
	a[N] = &a[0][0];
	a[N+1] = &a[0][m*o+1];
    
    //
    
	mpfr_t **atest = new mpfr_t * [N+2];
	atest[0] = new mpfr_t [(N+2)*(m*o+1)];
	
	for (int i=1; i<N+2; i++){
		atest[i] = &atest[0][i*(m*o+1)];
	}
	
	atest[N] = &atest[0][0];
	atest[N+1] = &atest[0][m*o+1];
    
    // creating arrays
    
    // variables
    
    // initialization
    
    mpfr_init(c);
    
    mpfr_init(zero);
    
    mpfr_init(par);
    
    mpfr_init2(w,1000);
    
    mpfr_init(actija);
    
    mpfr_init(actijt);
    
    mpfr_init(action);
    
    for (int k=0; k<o+1; k++){
		for (int l=0; l<o+1; l++){
			mpfr_init(A[k][l]);
        }
	}
    
    for (int k=0; k<o+1; k++){
		for (int l=0; l<o+1; l++){
			mpfr_init(M[k][l]);
		}
	}
    
    for (int p=0; p<N+2; p++){
		for (int l=0; l<m*o+1; l++){
			mpfr_init(a[p][l]);
		}
	}
    
    for (int p=0; p<N+2; p++){
		for (int l=0; l<m*o+1; l++){
			mpfr_init(atest[p][l]);
		}
	}
    
    // initialization
    
    // calculating the zero
    
    mpfr_set_flt(zero,1.0,MPFR_RNDD);
    
    char * ch=new char [1];
    
    mpfr_exp_t exp;
    
	int expzero=145;
	
    mpfr_set_flt(zero,10.0,MPFR_RNDD);
    
    mpfr_pow_si(zero,zero,-expzero,MPFR_RNDD);
    
    // calculating the zero
    
    // calculating the A coefficients
    
    mpfr_set_ui(c,1,MPFR_RNDD);
    mpfr_div_ui(c,c,3,MPFR_RNDD);
    
    mpfr_set_ui(A[0][0],2,MPFR_RNDD);
    mpfr_add(A[0][0],A[0][0],c,MPFR_RNDD);
    
    mpfr_set_ui(A[1][1],5,MPFR_RNDD);
    mpfr_add(A[1][1],A[1][1],c,MPFR_RNDD);
    
    mpfr_set(A[2][2],A[0][0],MPFR_RNDD);
    
    mpfr_set_si(A[0][1],-3,MPFR_RNDD);
    
    mpfr_add(A[0][1],A[0][1],c,MPFR_RNDD);
    
    mpfr_set(A[1][0],A[0][1],MPFR_RNDD);
    mpfr_set(A[2][1],A[0][1],MPFR_RNDD);
    mpfr_set(A[1][2],A[0][1],MPFR_RNDD);
    
    mpfr_set(A[0][2],c,MPFR_RNDD);
    mpfr_set(A[2][0],c,MPFR_RNDD);
    
    // calculating the A coefficients
    
    // calculating the M coefficients
    
    mpfr_set_d(M[0][0],1.0,MPFR_RNDD);
    mpfr_set_d(M[0][1],0.0,MPFR_RNDD);
    mpfr_set_d(M[0][2],0.0,MPFR_RNDD);
    mpfr_set_d(M[1][0],-3.0,MPFR_RNDD);
    mpfr_set_d(M[1][1],4.0,MPFR_RNDD);
    mpfr_set_d(M[1][2],-1.0,MPFR_RNDD);
    mpfr_set_d(M[2][0],2.0,MPFR_RNDD);
    mpfr_set_d(M[2][1],-4.0,MPFR_RNDD);
    mpfr_set_d(M[2][2],2.0,MPFR_RNDD);
    
    // calculating the M coefficients
    
    // saturating a to width w
    
    set_init(a,N,k,m,o,J,prec);
    set_init(atest,N,k,m,o,J,prec);
	
    
    int i=1;
    int j=1;
    int numrej=0;
    int plus;
    int minus;
    int count;
    
    
    mpfr_set_flt(w,50.0,MPFR_RNDU);
    mpfr_pow_si(w,w,-1,MPFR_RNDU);
    
//    double w=0.0001;

///// insert here from while (numrej<numvar) {
	
	while (numrej<numvar) {
        
        plus=0;
        minus=0;
        count=0;
        
        while (plus==0) {
            
            mpfr_add(atest[i][j],a[i][j],w,MPFR_RNDD);
			            
            action2ij(a,&actija, N, m, o, i, j, J, A, M, zero, prec);
            action2ij(atest,&actijt, N, m, o, i, j, J, A, M, zero, prec);
            			
            if (mpfr_cmp(actijt,actija)<0) {
                mpfr_set(a[i][j],atest[i][j],MPFR_RNDD);
                count=count+1;
            }
            else
            {
                mpfr_set(atest[i][j],a[i][j],MPFR_RNDD);
                plus=1;
            }
        }
        if (count==0)
        {
            while (minus==0) {
                
                mpfr_sub(atest[i][j],a[i][j],w,MPFR_RNDD);
                action2ij(a,&actija, N, m, o, i, j, J, A, M, zero, prec);
                action2ij(atest,&actijt, N, m, o, i, j, J, A, M, zero, prec);
                				
                
                if (mpfr_cmp(actijt,actija)<0) {
                    mpfr_set(a[i][j],atest[i][j],MPFR_RNDD);
                    count=count+1;
					
                }
                else
                {
                    mpfr_set(atest[i][j],a[i][j],MPFR_RNDD);
                    minus=1;
                }
            }
        }
        if (count==0) {
            numrej=numrej+1;
        }
        else
        {
            numrej=0;
        }
        if (j<maxj) {
            j=j+1;
        }
        else
        {
            j=1;
            if (i<N) {
                i=i+1;
            }
            else
            {
                i=1;
            }
        }
    }
    
///// insert here from while (numrej<numvar) }
    
    putchar('\n');
    
    action2(a,&action, N, m, o, J, A, M, zero, prec);
    
    mpfr_out_str (stdout, 10, 0, action, MPFR_RNDD);
    
    putchar('\n');
    
    
    i=3;
    
    j=9;
    
    putchar('\n');

    
    // saturating a to width w
    

    delete[] ch;
    delete[] A[0];
    delete[] A;
    delete[] M[0];
    delete[] M;
    delete[] a[0];
    delete[] a;
    delete[] atest[0];
    delete[] atest;
    
    mpfr_clear(zero);
    mpfr_clear(c);
    mpfr_clear(w);
    mpfr_clear(actija);
    mpfr_clear(actijt);
    
    return 0;
}
