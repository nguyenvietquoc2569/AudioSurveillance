/*
 * Panda_GMM.h
 *
 *  Created on: Jan 1, 2013
 *      Author: pandaubu
 */

#ifndef PANDA_GMM_H_
#define PANDA_GMM_H_




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>

typedef struct Panda_GMM_Component
{
	double N;
	double pb;
	double* Mu;
	double** R;
    double** invR;
	double Const;

} Panda_GMM_Component;

typedef struct Panda_GMM
{
	int NumOfMix;
	int dim;
	double Rmin;
	Panda_GMM_Component *Coms;
	double rissanen;
	double loglikelihood;
	double **pnk;
	int pns;
    char name[100];
    double likelihood;

} Panda_GMM ;

typedef struct Panda_GMM_Data
{
	double ** data;
	int dim;
	int nSamples;
} Panda_GMM_Data;

Panda_GMM *Panda_GMM_Init(Panda_GMM_Data *data,int initK);
void Panda_GMM_Normalize(Panda_GMM *GMM);
void Panda_GMM_Normalize_withoutPb(Panda_GMM *GMM);
void Panda_GMM_EMIterate(Panda_GMM_Data  *data,Panda_GMM *GMM);
double Panda_GMM_EStep(Panda_GMM_Data  *data,Panda_GMM *GMM);
void Panda_GMM_MStep(Panda_GMM_Data  *data,Panda_GMM *GMM);
void Panda_GMM_FreeGMM(Panda_GMM *GMM);
void Panda_GMM_FreeComponent(Panda_GMM_Component *Com,int dim);
Panda_GMM *Panda_GMM_Util_DuplicateGMM(Panda_GMM *GMM);
Panda_GMM_Component *Panda_GMM_Util_DuplicateGMMComponent(Panda_GMM_Component *Com,int dim);
Panda_GMM_Component *Panda_GMM_Component_Combine(Panda_GMM_Component *Com1,Panda_GMM_Component *Com2,int dim);
double Panda_GMM_Component_Distance(Panda_GMM_Component *Com1,Panda_GMM_Component *Com2,int dim);
Panda_GMM *Panda_GMM_MDLReduceOrder(Panda_GMM *GMM);
Panda_GMM *Panda_GMM_Create(Panda_GMM_Data *data,int initK,int finalK);
double Panda_GMM_likelihood(Panda_GMM_Data *data,Panda_GMM *GMM);
double *Panda_GMM_Gaupdf(Panda_GMM_Data *data,Panda_GMM *GMM);
///Update
Panda_GMM *Panda_GMM_IncreamentalLearning(Panda_GMM *GMM,Panda_GMM_Data *data);
double* Panda_GMM_IncreamentalLearning_Pix(Panda_GMM *GMM,Panda_GMM_Data *data);
void Panda_GMM_IncreamentalLearning_PlusPnk(Panda_GMM *GMM,double **pnk,int nsam);
double* Panda_GMM_IncreamentalLearning_SumPnk(double **pnk,int nsam,int NumMix);
double** Panda_GMM_IncreamentalLearning_pnk(Panda_GMM *GMM,Panda_GMM_Data *data);
//save file
void Panda_GMM_Save(char *filepath,Panda_GMM* GMM);
Panda_GMM* Panda_GMM_Read(char *filepath);
Panda_GMM* Panda_GMM_Load_backgroundModel(char *filepath);
#endif /* PANDA_GMM_H_ */
