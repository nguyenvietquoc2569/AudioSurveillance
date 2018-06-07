#include <stdlib.h>
#include <math.h>

#define PI  3.1415926535897932

double* DCT_II(double* x, int N)
{
	double* X = (double*) malloc(sizeof(double) * N);
	int k,n;
	
	for(k = 0; k < N; k++) 
	{
        X[k] = 0.0;
        
        for(n = 0; n < N; n++) 
        {
            X[k] = (X[k] + (x[n]*cos((PI/N)*(n-.5)*k)));
        }
    }
    
    return X;
}
