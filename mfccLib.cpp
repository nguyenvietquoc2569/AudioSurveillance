#define PI M_PI
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "fft.h"
#include "WavTool.h"
#include "mfccLib.h"
#include "utils.h"

//#include "dct.h"
//#include "complex_simple.h"
/*
 *
 */
double Panda_support_M2F(double Mel)
{
    return 700*exp(Mel/1127)-700;
}
double Panda_support_F2M(double Hz)
{
    return 1127*log(1+Hz/700);
}

double* Panda_MfccParameter_createWindow(int length,char* name)
{
    double* wi;
    int i;
    wi = (double*) malloc(sizeof(double)*length);

        if (strcmp(name,"hann")==0)
        {
            for(i=0;i<length;i++)
            {

                wi[i]=0.5*(1-cos(2*M_PI*(i)/(length-1)));
                //printf("%f : ",wi[i]);
            }
        }
        else
        {
            for(i=0;i<length;i++)
            {

                wi[i]=0.54f - 0.46*cos(2*M_PI*((double)i)/((double)(length-1)));

            }
        }

    return wi;
}
double** Panda_MfccParameter_createDCTmatrix(int N,int M)
{
    double** wi;
    double* te;
    double sq;
    int i,j;

    sq =sqrt((double)2/M);
    wi = (double**)malloc(sizeof(double*)*N);
    te= (double*)malloc(sizeof(double)*M);
    for(i=1;i<=M;i++)
    {
        te[i-1]=M_PI*(i-0.5)/M;
        //printf("%f \n",te[i-1]);
    }
    for(i=0;i<N;i++)
    {
        wi[i]= (double*)malloc(sizeof(double)*M);
        for(j=0;j<M;j++)
        {
            wi[i][j]=sq*cos(i*te[j]);
            /*printf("%f  ",wi[i][j]);
            printf("%d:%d:%f   ...... ",i,j,wi[i][j]);*/
        }
    }
    free(te);
    return wi;
}
double** Panda_MfccParameter_createTribandmatrix(int M,int K,long f_low,long f_high,long fs)
{

    double **wi;
    int i,j;
    double f_min,f_max,temp,*f,*fw,*c,*cw;

    f=(double*)malloc(sizeof(double)*K);
    fw=(double*)malloc(sizeof(double)*K);
    f_min=0;f_max=0.5f * fs;
    temp=(double)(f_max-f_min)/(K-1);

    for(i=0;i<K;i++)
    {
        f[i]=f_min+i*temp;
        fw[i]=Panda_support_F2M(f[i]);
    }

    c=(double*)malloc(sizeof(double)*(M+2));
    cw=(double*)malloc(sizeof(double)*(M+2));
    temp=(Panda_support_F2M(f_high)-Panda_support_F2M(f_low))/(M+1);
    for (i=0;i<=M+1;i++)
    {
        cw[i]=Panda_support_F2M(f_low)+i*temp;
        c[i]=Panda_support_M2F(cw[i]);
    }

    /*for(i=0;i<M+3;i++) //kiem tra
    {
        printf("%d:%f   ",i,c[i]);
    }*/

    wi = (double**)malloc(sizeof(double*)*M);
    for(i=0;i<M;i++)
    {
        wi[i]=(double*)malloc(sizeof(double)*K);
        memset(wi[i],0,sizeof(double)*K);
        for(j=0;j<K;j++)
        {
            if (f[j]>=c[i] && f[j]<=c[i+1])
            {
                wi[i][j]=(f[j]-c[i])/(c[i+1]-c[i]);
            }
            if (f[j]>=c[i+1] && f[j]<=c[i+2])
            {
                wi[i][j]=(c[i+2]-f[j])/(c[i+2]-c[i+1]);
            }
        }

    }
    free(f);
    free(fw);
    free(c);
    free(cw);

    return wi;
}
double* PandaSampleToMfcc(Panda_mfccParameter *Param,short *speech,long NSpeech)
{

    int i,j;
    double sum;
    complex *sound_buffer;
    complex *fft;
    double *mag;
    double* FBE;
    double *CC;
    if (Param->NSframe != NSpeech)
    {
        printf("size frames is not the same with parameter");
        return NULL;
    }




    sound_buffer = (complex *)malloc(sizeof(complex)*Param->nfft);


    for(i=1;i<NSpeech;i++)
    {
    	sound_buffer[i].re=((double)(speech[i] -  speech[i-1]*Param->alpha))*Param->window[i];
        //printf("%d:%f\n ",i+1,sound_buffer[i].re);
        sound_buffer[i].im=0;
    }
    sound_buffer[0].re=((double)(speech[0]))*Param->window[0];
    sound_buffer[0].im=0;

    fft = FFT_simple(sound_buffer,Param->nfft);
    free(sound_buffer);

    mag = (double*)malloc(sizeof(double)*(Param->nfft/2+1));
    for(i=0;i<(Param->nfft/2+1);i++)
    {
        mag[i]=complex_magnitude(fft[i]);
        //printf("%d:%f\n ",i+1,fft[i].re);
    }
    free(fft);



    FBE= (double *)malloc(sizeof(double)*Param->M);
    for(j=0;j<Param->M;j++)
    {
        sum =0;
        for(i=0;i<(Param->nfft/2+1);i++)
        {
            sum=sum+Param->tribank[j][i]*mag[i];
        }
        FBE[j]=log(sum);
    }
    free(mag);

    CC= (double *)malloc(sizeof(double)*Param->N);
    for(j=0;j<Param->N;j++)
    {
        sum =0;
        for(i=0;i<=Param->M;i++)
        {
            sum=sum+Param->DCT[j][i]*FBE[i];
        }
        CC[j]=sum;
        /* lift mfcc */
        CC[j]=CC[j]*Param->ceplifter[j];

    }
    free(FBE);






    return CC;
}
double* PandaSampleToMfccWithoutA(Panda_mfccParameter *Param,double *speech,long NSpeech)
{

    int i,j;
    double sum;
    complex *sound_buffer;
    complex *fft;
    double *mag;
    double* FBE;
    double *CC;
    if (Param->NSframe != NSpeech)
    {
        printf("size frames is not the same with parameter");
        return NULL;
    }




    sound_buffer = (complex *)malloc(sizeof(complex)*Param->nfft);
    for(i=1;i<Param->nfft;i++)
    {
    	sound_buffer[i].re=0;
    	sound_buffer[i].im=0;
    }

    for(i=1;i<NSpeech;i++)
    {
    	sound_buffer[i].re=((double)(speech[i]))*Param->window[i];
        //printf("%d:%f\n ",i+1,sound_buffer[i].re);
        sound_buffer[i].im=0;
    }
    sound_buffer[0].re=((double)(speech[0]))*Param->window[0];
    sound_buffer[0].im=0;

    fft = FFT_simple(sound_buffer,Param->nfft);
    free(sound_buffer);

    mag = (double*)malloc(sizeof(double)*(Param->nfft/2+1));
    for(i=0;i<(Param->nfft/2+1);i++)
    {
        mag[i]=complex_magnitude(fft[i]);
        //printf("%d:%f\n ",i+1,fft[i].re);
    }
    free(fft);



    FBE= (double *)malloc(sizeof(double)*Param->M);
    for(j=0;j<Param->M;j++)
    {
        sum =0;
        for(i=0;i<(Param->nfft/2+1);i++)
        {
            sum=sum+Param->tribank[j][i]*mag[i];
        }
        FBE[j]=log(sum);
    }
    free(mag);

    CC= (double *)malloc(sizeof(double)*Param->N);
    for(j=0;j<Param->N;j++)
    {
        sum =0;
        for(i=0;i<=Param->M;i++)
        {
            sum=sum+Param->DCT[j][i]*FBE[i];
        }
        CC[j]=sum;
        /* lift mfcc */
        CC[j]=CC[j]*Param->ceplifter[j];

    }
    free(FBE);






    return CC;
}
Panda_mfccParameter* Panda_Mfcc_CreateParam(long fs,long NS_frame,long NS_ShiftFrame ,char* window,long f_low,long f_high,int M,int N,int L,double alpha)
{
	int i;
    Panda_mfccParameter* para;
    para = (Panda_mfccParameter*)malloc(sizeof(Panda_mfccParameter));

    para->samplerate =fs;
    para->NSframe =NS_frame;
    para->nfft= pow(2,round(log2(NS_frame)+0.4999f));
    para->M = M;
    para->N = N;
    para->alpha=alpha;

    para->tribank = Panda_MfccParameter_createTribandmatrix(M,(para->nfft/2)+1,f_low,f_high,fs);
    para->window=Panda_MfccParameter_createWindow(para->NSframe,window);
    para->DCT=Panda_MfccParameter_createDCTmatrix(N,M);

    para->ceplifter = (double*)malloc(sizeof(double)*N);
    for(i=0;i<N;i++)
    {
    	para->ceplifter[i]=(1+0.5*L*sin(M_PI*((double)i)/L));
    }
    para->NS_ShiftFrame=NS_ShiftFrame;

    return para;
}


int Panda_Mfcc_readFileWav(double*** CC,Panda_mfccParameter *q, char * filepath)
{
	wav_inf* file;
	double **re,*temp;
	file = Wavtool_readSam(filepath);
	//Wavtool_printinfo(filepath);
	int le,i;
	long s;


	le = (file->nsample-q->NSframe)/q->NS_ShiftFrame + 1;
	re = (double**)malloc(sizeof(double*)*le);
	
	temp=(double*)malloc(sizeof(double)*file->nsample);
	for(i=file->nsample-1;i>=1;i--)
    {
        temp[i]=((double)(file->samples[i])) - ((double)file->samples[i-1])*q->alpha;
        //printf("%d:%f\n ",i+1,sound_buffer[i].re);
    }
    temp[0]=(double)file->samples[0];
        //printf("\nQ1 %d....................\n",le);
    s=0;
    for(i=0;i<le;i++)
    {
        re[i]=PandaSampleToMfccWithoutA(q,&(temp[s]),q->NSframe);
        //printf("\n %f \n",re[i][0]);
        s=s+q->NS_ShiftFrame;
        //printf("\nQ1 %d....................\n",i);
    }


    free(temp);
    free(file);
    *CC = re;

	return le;
}


int Panda_Mfcc_readForder (double*** CC,Panda_mfccParameter *q, char * path)
{
	int nmf,i;
	nmf=0;
	double **C,*temp;
	char fullstr[200];
    // first off, we need to create a pointer to a directory
    DIR *pdir = NULL; // remember, it's good practice to initialise a pointer to NULL!

    wav_hdr *hrd;
    wav_inf *file;
    long nsample;
    long dem,s;
    dem=-1;


    pdir = opendir (path); // "." will refer to the current directory
    struct dirent *pent = NULL;
    if (pdir == NULL) // if pdir wasn't initialised correctly
    { // print an error message and exit the program
        printf ("\nERROR! pdir could not be initialised correctly");
        return nmf; // exit the function
    } // end if

    while (pent = readdir (pdir)) // while there is still something in the directory to list
    {
        if (pent == NULL) // if pent has not been initialised correctly
        { // print an error message, and exit the program
            printf ("\nERROR! pent could not be initialised correctly");
            return nmf; // exit the function
        }
        if (strcmp(pent->d_name,"..")==0 || strcmp(pent->d_name,".")==0) continue;
        strcpy(fullstr,path);
        strcat(fullstr,pent->d_name);
        hrd = Wavtool_read(fullstr);
        nsample = hrd->Subchunk2Size/(hrd->bitsPerSample/8);
        nmf = nmf + (nsample-q->NSframe)/q->NS_ShiftFrame + 1;
        // otherwise, it was initialised correctly. let's print it on the console:
        printf ("%s\n", pent->d_name);
    }

    // finally, let's close the directory
    closedir (pdir);

    C = (double**)malloc(sizeof(double*)*nmf);


    pdir = opendir (path);
    while (pent = readdir (pdir)) // while there is still something in the directory to list
        {
            if (pent == NULL) // if pent has not been initialised correctly
            { // print an error message, and exit the program
                printf ("\nERROR! pent could not be initialised correctly");
                return nmf; // exit the function
            }
            if (strcmp(pent->d_name,"..")==0 || strcmp(pent->d_name,".")==0) continue;
            strcpy(fullstr,path);
            strcat(fullstr,pent->d_name);

            file = Wavtool_readSam(fullstr);
            temp=(double*)malloc(sizeof(double)*file->nsample);
            for(i=file->nsample-1;i>=1;i--)
            {
            	temp[i]=((double)(((double)file->samples[i])- ((double)file->samples[i-1])*q->alpha));
                //printf("%d:%f\n ",i+1,sound_buffer[i].re);
            }
            temp[0]=(double)file->samples[0];
            s=0;
            while (s+q->NSframe-1<file->nsample)
            {
            	dem++;
            	C[dem]=PandaSampleToMfccWithoutA(q,&(temp[s]),q->NSframe);


            	//printf("\n %f \n",re[i][0]);
            	s=s+q->NS_ShiftFrame;
            }
            free(temp);
            free(file);

            // otherwise, it was initialised correctly. let's print it on the console:
            printf ("%s\n", pent->d_name);
    }
	*CC=C;
    return nmf;
}
