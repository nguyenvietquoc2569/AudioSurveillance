/*
 * Panda_GMM_Util.h
 *
 *  Created on: Jan 1, 2013
 *      Author: pandaubu
 */

#ifndef PANDA_GMM_UTIL_H_
#define PANDA_GMM_UTIL_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include <float.h>


double Panda_GMM_Util_Cov_Var(double *mean,double **x,int dim,int num,int posx,int posy);
double Panda_GMM_Util_Cov_Var_NorWithN(double *mean,double **x,int dim,int num,int posx,int posy);
double **Panda_GMM_Util_Cov_NormalWithN(double **x,int n,int dim);
double **Panda_GMM_Util_Cov(double **x,int n,int dim);
double *Panda_GMM_Util_Diag(double **x,int n);
void Panda_GMM_Util_freedouble(double ***x,int n);
double Panda_GMM_Util_Mean(double *x,int n);
double **Panda_GMM_Util_Eye(int n);
double **Panda_GMM_Util_Copy2(double **x,int n);
double **Panda_GMM_Util_Copy3(double **x,int n,int m);
double *Panda_GMM_Util_Copy(double *x,int n);
double Panda_GMM_Util_Determinant(double **a,int n);
double **Panda_GMM_Util_Zeros(int n,int m);
double **Panda_GMM_Util_Ones(int n,int m);
double *Panda_GMM_Util_Multi_vecAma(double *vec,double** mat,int n_m,int m_m,double scl);
double Panda_GMM_Util_Multi_vecAvec(double *vec1,double *vec2,int n);
double *Panda_GMM_Util_Minus_vecAvec(double *vec1,double *vec2,int n);

#endif /* PANDA_GMM_UTIL_H_ */

