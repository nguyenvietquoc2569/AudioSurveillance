/*
 * mfccLib.h
 *
 *  Created on: Jan 11, 2013
 *      Author: panda
 */

#ifndef MFCCLIB_H_
#define MFCCLIB_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "fft.h"
#include "WavTool.h"
//#include "dct.h"
//#include "complex_simple.h"
/*
 *
 */
double Panda_support_M2F(double Mel);
double Panda_support_F2M(double Hz);
typedef struct Panda_MfccParameter
{
    long samplerate;
    int nfft;
    double* window;
    long NSframe;
    long NS_ShiftFrame;
    double **tribank;
    double **DCT;
    int M,N;
    double *ceplifter;
    double alpha;


} Panda_mfccParameter;

double* Panda_MfccParameter_createWindow(int length,char* name);
double** Panda_MfccParameter_createDCTmatrix(int N,int M);
double** Panda_MfccParameter_createTribandmatrix(int M,int K,long f_low,long f_high,long fs);
double* PandaSampleToMfcc(Panda_mfccParameter *Param,short *speech,long NSpeech);
Panda_mfccParameter* Panda_Mfcc_CreateParam(long fs,long NS_frame,long NS_ShiftFrame ,char* window,long f_low,long f_high,int M,int N,int L,double alpha);
int main1();
int Panda_Mfcc_readFileWav(double*** CC,Panda_mfccParameter *q, char * filepath);
double* PandaSampleToMfccWithoutA(Panda_mfccParameter *Param,double *speech,long NSpeech);
int Panda_Mfcc_readForder (double*** CC,Panda_mfccParameter *q, char * path);
#endif /* MFCCLIB_H_ */
