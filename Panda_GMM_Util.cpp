#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include <float.h>


double Panda_GMM_Util_Cov_Var(double *mean,double **x,int dim,int num,int posx,int posy)
{
	double re;
	int i;
	re =0 ;
	for (i=0;i<num;i++)
	{
		re = re + (mean[posx]-x[i][posx])*(mean[posy]-x[i][posy]);
	}
	re=re/((double)num-1.0f);
	return re;
}
double Panda_GMM_Util_Cov_Var_NorWithN(double *mean,double **x,int dim,int num,int posx,int posy)
{
	double re;
	int i;
	re =0 ;
	for (i=0;i<num;i++)
	{
		re = re + (mean[posx]-x[i][posx])*(mean[posy]-x[i][posy]);
	}
	re=re/((double)num);
	return re;
}
double **Panda_GMM_Util_Cov(double **x,int n,int dim)
{
	double **re,*mean;
	int i,j;


	re = (double**)malloc(sizeof(double*)*dim);
	mean = (double*)malloc(sizeof(double)*dim);

	for(i=0;i<dim;i++)
	{
		re[i] = (double*)malloc(sizeof(double)*dim);
		mean[i]=0;
	}

	//Panda calcs Mean Vector
	for(i=0;i<n;i++)
	{
		for(j=0;j<dim;j++)
		{
			mean[j]=mean[j]+x[i][j];
		}
	}

	for(i=0;i<dim;i++)
		mean[i]=mean[i]/n;

	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
			{
				re[j][i]=Panda_GMM_Util_Cov_Var(mean,x,dim,n,i,j);
				re[i][j]=re[j][i];

			}
	}
	free(mean);
	return re;
}
double **Panda_GMM_Util_Cov_NormalWithN(double **x,int n,int dim)
{
	double **re,*mean;
	int i,j;


	re = (double**)malloc(sizeof(double*)*dim);
	mean = (double*)malloc(sizeof(double)*dim);

	for(i=0;i<dim;i++)
	{
		re[i] = (double*)malloc(sizeof(double)*dim);
		mean[i]=0;
	}

	//Panda calcs Mean Vector
	for(i=0;i<n;i++)
	{
		for(j=0;j<dim;j++)
		{
			mean[j]=mean[j]+x[i][j];
		}
	}

	for(i=0;i<dim;i++)
        mean[i]=mean[i]/(double)n;

	for(i=0;i<dim;i++)
	{
		for(j=i;j<dim;j++)
			{
				re[j][i]=Panda_GMM_Util_Cov_Var_NorWithN(mean,x,dim,n,i,j);
				re[i][j]=re[j][i];
			}
	}
	free(mean);
	return re;
}
double *Panda_GMM_Util_Diag(double **x,int n)
{
	double *re;
	int i;
	re = (double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++) re[i]=x[i][i];
	return re;
}
void Panda_GMM_Util_freedouble(double ***x,int n)
{
	int i;
	for(i=0;i<n;i++) free((*x)[i]);
	free(*x);
}
double Panda_GMM_Util_Mean(double *x,int n)
{
	double re;
	int i;
	re=0;
	for (i=0;i<n;i++)
	{
		re=re+x[i];
	}
	re=(double)re/(double)n;
	return re;
}
double **Panda_GMM_Util_Eye(int n)
{
	double **re;
	int i,j;

	re= (double**)malloc(sizeof(double*)*n);
	for (i=0;i<n;i++)
	{
		re[i]=(double*)malloc(sizeof(double)*n);
	}
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++) re[i][j]=0;
	}
	for (i=0;i<n;i++)
		{
			re[i][i]=1;
		}


	return re;
}
double *Panda_GMM_Util_Copy(double *x,int n)
{
	double *re;
	int i;
	re=(double*)malloc(sizeof(double)*n);
	for (i=0;i<n;i++) re[i]=x[i];
	return re;
}
double **Panda_GMM_Util_Copy2(double **x,int n)
{
	double **re;

	int i,j;

	re=(double**)malloc(sizeof(double*)*n);
	for (i=0;i<n;i++) re[i]=(double*)malloc(sizeof(double)*n);

	for (i=0;i<n;i++) for (j=0;j<n;j++) re[i][j]=x[i][j];
	return re;
}
double **Panda_GMM_Util_Copy3(double **x,int n,int m)
{
	double **re;

	int i,j;

	re=(double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
    {
        re[i]=(double*)malloc(sizeof(double)*m);
        memcpy(re[i],x[i],m*sizeof(double));
    }


	return re;
}
double Panda_GMM_Util_Determinant(double **a,int n)
{
    int i,j,j1,j2 ;                    // general loop and matrix subscripts
    double det = 0 ;                   // init determinant
    double **m = NULL ;                // pointer to pointers to implement 2d
                                       // square array

    if (n < 1)    {   }                // error condition, should never get here

    else if (n == 1) {                 // should not get here
        det = a[0][0] ;
        }

    else if (n == 2)  {                // basic 2X2 sub-matrix determinate
                                       // definition. When n==2, this ends the
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1] ;// the recursion series
        }


                                       // recursion continues, solve next sub-matrix
    else {                             // solve the next minor by building a
                                       // sub matrix
        det = 0 ;                      // initialize determinant of sub-matrix

                                           // for each column in sub-matrix
        for (j1 = 0 ; j1 < n ; j1++) {
                                           // get space for the pointer list
            m = (double **) malloc((n-1)* sizeof(double *)) ;

            for (i = 0 ; i < n-1 ; i++)
                m[i] = (double *) malloc((n-1)* sizeof(double)) ;
                       //     i[0][1][2][3]  first malloc
                       //  m -> +  +  +  +   space for 4 pointers
                       //       |  |  |  |          j  second malloc
                       //       |  |  |  +-> _ _ _ [0] pointers to
                       //       |  |  +----> _ _ _ [1] and memory for
                       //       |  +-------> _ a _ [2] 4 doubles
                       //       +----------> _ _ _ [3]
                       //
                       //                   a[1][2]
                      // build sub-matrix with minor elements excluded
            for (i = 1 ; i < n ; i++) {
                j2 = 0 ;               // start at first sum-matrix column position
                                       // loop to copy source matrix less one column
                for (j = 0 ; j < n ; j++) {
                    if (j == j1) continue ; // don't copy the minor column element

                    m[i-1][j2] = a[i][j] ;  // copy source element into new sub-matrix
                                            // i-1 because new sub-matrix is one row
                                            // (and column) smaller with excluded minors
                    j2++ ;                  // move to next sub-matrix column position
                    }
                }

            det += pow(-1.0,1.0 + j1 + 1.0) * a[0][j1] * Panda_GMM_Util_Determinant(m,n-1) ;
                                            // sum x raised to y power
                                            // recursively get determinant of next
                                            // sub-matrix which is now one
                                            // row & column smaller

            for (i = 0 ; i < n-1 ; i++) free(m[i]) ;// free the storage allocated to
                                            // to this minor's set of pointers
            free(m) ;                       // free the storage for the original
                                            // pointer to pointer
        }
    }
    return det ;
}
double **Panda_GMM_Util_Zeros(int n,int m)
{
	double **re;
	int i,j;
	re=(double**)malloc(sizeof(double*)*n);
	for(i=0;i<n;i++)
	{
		re[i]=(double *)malloc(sizeof(double)*m);
        memset(re[i],0,sizeof(double)*m);
	}

	return re;
}
double **Panda_GMM_Util_Ones(int n,int m)
{
	double **re;
	int i,j;
	re=(double**)malloc(sizeof(double*)*n);
	for(i=0;i<n;i++)
	{
		re[i]=(double *)malloc(sizeof(double)*m);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
			re[i][j]=1;
	}

	return re;
}
double *Panda_GMM_Util_Multi_vecAma(double *vec,double** mat,int n_m,int m_m,double scl) //vec have size like n_m
{
	double* re;
	int i,j ;
	re = (double*)malloc(sizeof(double)*m_m);
	for(i=0;i<m_m;i++)
	{
		re[i]=0;
		for(j=0;j<n_m;j++)
			re[i]=re[i]+mat[j][i]*vec[j];
		re[i]=re[i]*scl;
	}
	return re;
}
double Panda_GMM_Util_Multi_vecAvec(double *vec1,double *vec2,int n)
{
	double re;
	int i;
	re = 0;
	for(i=0;i<n;i++)
	{
		re=re+vec1[i]*vec2[i];
	}

	return re;
}
double *Panda_GMM_Util_Minus_vecAvec(double *vec1,double *vec2,int n)
{
	double* re;
	int i ;
	re = (double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++)
	{
		re[i]=vec1[i]-vec2[i];
	}

	return re;
}



