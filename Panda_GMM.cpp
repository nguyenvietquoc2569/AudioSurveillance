/*
 * Panda_GMM.c
 *
 *  Created on: Dec 31, 2012
 *      Author: pandaubu
 */

#include "Panda_GMM.h"
#include "Panda_GMM_Util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include <string.h>
#include "invert.h"
#include "QDebug"

double matrix_det ( double **in_matrix, int n )
{
  int i, j, k;
  double **matrix;
  double det = 1,te;

  matrix = (double **)malloc(sizeof(double*)*n);

  for ( i = 0; i < n; i++ )
    matrix[i] = (double *)malloc(sizeof(double)*n);

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      matrix[i][j] = in_matrix[i][j];
  }

  for ( k = 0; k < n; k++ ) {
    if ( matrix[k][k] == 0 ) {
      int ok = 0;

      for ( j = k; j < n; j++ ) {
        if ( matrix[j][k] != 0 )
          ok = 1;
      }

      if ( ok==0 )
        return 0;

      for ( i = k; i < n; i++ )
      {
        te= matrix[i][j];matrix[i][j]= matrix[i][k] ;matrix[i][k]=te;
      }

      det = -det;
    }

    det *= matrix[k][k];

    if ( k + 1 < n ) {
      for ( i = k + 1; i < n; i++ ) {
        for ( j = k + 1; j < n; j++ )
          matrix[i][j] = matrix[i][j] - matrix[i][k] *
          matrix[k][j] / matrix[k][k];
      }
    }
  }

  for ( i = 0; i < n; i++ )
    free(matrix[i]);

  free(matrix);

  return det;
}

void Panda_GMM_Normalize(Panda_GMM *GMM)
{
	double sum;
	int i;
	sum =0;
	for(i=0;i<GMM->NumOfMix;i++)
	{
		sum = sum + GMM->Coms[i].pb;
	}

	for(i=0;i<GMM->NumOfMix;i++)
	{
		GMM->Coms[i].pb = GMM->Coms[i].pb/sum;

		Panda_GMM_Util_freedouble(&(GMM->Coms[i].invR),GMM->dim);
		GMM->Coms[i].invR = Panda_GMM_Util_Copy2(GMM->Coms[i].R,GMM->dim);
		invert(GMM->Coms[i].invR,GMM->dim);
        GMM->Coms[i].Const = -(((double)GMM->dim) * log(2*M_PI) +
                               ((double)log(matrix_det(GMM->Coms[i].R,GMM->dim))))/2.0f;
	}

}
void Panda_GMM_Normalize_withoutPb(Panda_GMM *GMM)
{
	int i;

	for(i=0;i<GMM->NumOfMix;i++)
	{


		Panda_GMM_Util_freedouble(&(GMM->Coms[i].invR),GMM->dim);
		GMM->Coms[i].invR = Panda_GMM_Util_Copy2(GMM->Coms[i].R,GMM->dim);
		invert(GMM->Coms[i].invR,GMM->dim);
		GMM->Coms[i].Const = -(((double)GMM->dim) * log(2*M_PI) + ((double)log(matrix_det(GMM->Coms[i].R,GMM->dim))))/2;//Panda_GMM_Util_Determinant
	}

}
Panda_GMM *Panda_GMM_Init(Panda_GMM_Data *data,int initK)
{
	Panda_GMM *re;
	double **R,**eye;
	double *diag;
	double period;
	int i,j,z,temp;


	re = (Panda_GMM*)malloc(sizeof(Panda_GMM));
	re->pns=0;
	eye=Panda_GMM_Util_Eye(data->dim);
	///////////////////////////xoa


	re->NumOfMix = initK;
	re->dim = data->dim;

    //sai ở hàm bên dưới :D
	R = Panda_GMM_Util_Cov_NormalWithN(data->data,data->nSamples,data->dim);

    //debug
    //Panda_debug::writeDouble2(R,data->dim,data->dim);


	diag = Panda_GMM_Util_Diag(R,data->dim);
	re->Rmin = Panda_GMM_Util_Mean(diag,data->dim)/((double)100000);
	if (re->Rmin==0)
	{
		re->Rmin=DBL_MIN;
	}

	free(diag);

	re->Coms = (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component)*initK);
//set com1
	re->Coms[0].N=0;
	re->Coms[0].pb = 1.0f/(double)initK;
	re->Coms[0].Mu = Panda_GMM_Util_Copy(data->data[0],data->dim);

	re->Coms[0].R = (double**)malloc(sizeof(double*)*data->dim);
	for(i=0;i<data->dim;i++)
		re->Coms[0].R[i] = (double*)malloc(sizeof(double)*data->dim);
	re->Coms[0].invR = (double**)malloc(sizeof(double*)*data->dim);
    for(i=0;i<data->dim;i++)
        re->Coms[0].invR[i] = (double*)malloc(sizeof(double)*data->dim);

	for(i=0;i<data->dim;i++)
		for(j=0;j<data->dim;j++)
			re->Coms[0].R[i][j]=R[i][j]+re->Rmin*eye[i][j];
//Set orther com
	period = ((double)data->nSamples-1)/((double)initK-1);
	for(z=1;z<initK;z++)
	{
		re->Coms[z].N=0;
		re->Coms[z].pb = 1.0f/(double)initK;
		temp=((int)roundf((z)*period+0.5)-1);
		re->Coms[z].Mu = Panda_GMM_Util_Copy(data->data[temp],data->dim);

		re->Coms[z].R = (double**)malloc(sizeof(double*)*data->dim);
		for(i=0;i<data->dim;i++)
			re->Coms[z].R[i] = (double*)malloc(sizeof(double)*data->dim);

		re->Coms[z].invR = (double**)malloc(sizeof(double*)*data->dim);
		for(i=0;i<data->dim;i++)
			re->Coms[z].invR[i] = (double*)malloc(sizeof(double)*data->dim);


		for(i=0;i<data->dim;i++)
			for(j=0;j<data->dim;j++)
					re->Coms[z].R[i][j]=R[i][j]+re->Rmin*eye[i][j];
	}


    //set pnk;
    re->pnk = (double**) malloc(sizeof(double*)*data->nSamples);
    for(i=0;i<data->nSamples;i++)
    {
        re->pnk[i] = (double*) malloc(sizeof(double)*initK);
    }

	Panda_GMM_Normalize(re);
	Panda_GMM_Util_freedouble(&R,data->dim);
	Panda_GMM_Util_freedouble(&eye,data->dim);

	return re;
}
void Panda_GMM_FreeGMM(Panda_GMM *GMM)
{
	int i;
    if (GMM->pnk != NULL)
        Panda_GMM_Util_freedouble(&(GMM->pnk),GMM->pns);
	GMM->pnk = NULL;
	for (i=GMM->NumOfMix-1;i>-1;i--)
	{
		Panda_GMM_Util_freedouble(&(GMM->Coms[i].R),GMM->dim);
		Panda_GMM_Util_freedouble(&(GMM->Coms[i].invR),GMM->dim);

		free(GMM->Coms[i].Mu);
	}

	free(GMM->Coms);
	free(GMM);
}
void Panda_GMM_FreeComponent(Panda_GMM_Component *Com,int dim)
{
	Panda_GMM_Util_freedouble(&(Com->R),dim);
	Panda_GMM_Util_freedouble(&(Com->invR),dim);
	free(Com->Mu);
    //free(Com);
}
Panda_GMM *Panda_GMM_Util_DuplicateGMM(Panda_GMM *GMM)
{
	int i,j,k;
	Panda_GMM *re;
	re= (Panda_GMM*) malloc(sizeof(Panda_GMM));
	re->rissanen = GMM->rissanen;
	re->loglikelihood = GMM->loglikelihood;
	re->dim = GMM->dim;
	re->Rmin = GMM->Rmin;
	re->NumOfMix = GMM->NumOfMix;
	re->pns = GMM->pns;
	re->Coms = (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component)*GMM->NumOfMix);
	for(i=0;i<GMM->NumOfMix;i++)
	{
		re->Coms[i].pb = GMM->Coms[i].pb;
		re->Coms[i].N = GMM->Coms[i].N;
		re->Coms[i].Const = GMM->Coms[i].Const;

		//copy R
		re->Coms[i].R = (double**)malloc(sizeof(double*)*GMM->dim);
		for(j=0;j<GMM->dim;j++)
		{
			re->Coms[i].R[j] = (double*)malloc(sizeof(double)*GMM->dim);
			for(k=0;k<GMM->dim;k++)
				re->Coms[i].R[j][k]=GMM->Coms[i].R[j][k];
		}
		//Copy invR
		re->Coms[i].invR = (double**)malloc(sizeof(double*)*GMM->dim);
		for(j=0;j<GMM->dim;j++)
		{
			re->Coms[i].invR[j] = (double*)malloc(sizeof(double)*GMM->dim);
			for(k=0;k<GMM->dim;k++)
				re->Coms[i].invR[j][k]=GMM->Coms[i].invR[j][k];
		}
		//Copy mean
		re->Coms[i].Mu = (double*)malloc(sizeof(double)*GMM->dim);
		for(j=0;j<GMM->dim;j++)
				re->Coms[i].Mu[j] = GMM->Coms[i].Mu[j];

		//Copy pnk
		if (GMM->pnk!=NULL)
		{
			re->pnk= Panda_GMM_Util_Copy3(GMM->pnk,GMM->pns,GMM->NumOfMix);
			re->pns = GMM->pns;
		}else
		{
			re->pns=0;
			re->pnk=NULL;
		}

	}
	return re;
}
Panda_GMM_Component *Panda_GMM_Util_DuplicateGMMComponent(Panda_GMM_Component *Com,int dim)
{
	Panda_GMM_Component *re;
	int j,k;
	re= (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component));
	re->pb = Com->pb;
	re->N = Com->N;
	re->Const = Com->Const;

	//copy R
	re->R = (double**)malloc(sizeof(double*)*dim);
	for(j=0;j<dim;j++)
	{
		re->R[j] = (double*)malloc(sizeof(double)*dim);
        memcpy(re->R[j],Com->R[j],sizeof(double)*dim);
//		for(k=0;k<dim;k++)
//			re->R[j][k]=Com->R[j][k];
	}
	//Copy invR
	re->invR = (double**)malloc(sizeof(double*)*dim);
	for(j=0;j<dim;j++)
	{
		re->invR[j] = (double*)malloc(sizeof(double)*dim);
        memcpy(re->invR[j],Com->invR[j],sizeof(double)*dim);
//		for(k=0;k<dim;k++)
//			re->invR[j][k]=Com->invR[j][k];
	}
	//copy Mean
	re->Mu = (double*)malloc(sizeof(double)*dim);
    memcpy(re->Mu,Com->Mu,sizeof(double)*dim);

//	for(j=0;j<dim;j++)
//		re->Mu[j] = Com->Mu[j];
	return re;
}
void Panda_GMM_Component_Copy(Panda_GMM_Component *To,Panda_GMM_Component *From,int dim)
{
	int i,j;
	To->Const=From->Const;
	To->N = From->N;
	To->pb = From->pb;
	for(i=0;i<dim;i++)
	{
		To->Mu[i]=From->Mu[i];
		for(j=0;j<dim;j++)
			{
				To->R[i][j]=From->R[i][j];
				To->invR[i][j]=From->invR[i][j];
			}
	}
}
void Panda_GMM_Component_Init(Panda_GMM_Component *C,int dim)
{
	int i;
	//C=(Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component));
	C->Mu = (double*)malloc(sizeof(double)*dim);
	C->R = (double**)malloc(sizeof(double*)*dim);
	C->invR = (double**)malloc(sizeof(double*)*dim);
	for(i=0;i<dim;i++)
	{
		C->R[i]=(double*)malloc(sizeof(double)*dim);
		C->invR[i]=(double*)malloc(sizeof(double)*dim);
	}
}
double Panda_GMM_Extend_Process1(double** pnk,int nsam,int ncom,Panda_GMM *GMM)
{
    double *ss,likelihood;
	int i,j;

	ss=(double*)malloc(sizeof(double)*nsam);



	for(i=0;i<nsam;i++)
	{
		ss[i]=0;


		for(j=0;j<ncom;j++)
		{
                    pnk[i][j]=((double)exp(((double)(pnk[i][j]))))*GMM->Coms[j].pb;
					ss[i]=ss[i]+pnk[i][j];
		}

	}

	likelihood = 0;

	//save pnk
	for(i=0;i<nsam;i++)
	{
        likelihood = likelihood + log(ss[i]);
		for(j=0;j<ncom;j++)
		{
			pnk[i][j]=pnk[i][j]/ss[i];
		}
	}
	if (GMM->pnk != NULL) Panda_GMM_Util_freedouble(&(GMM->pnk),nsam);
	GMM->pnk = Panda_GMM_Util_Copy3(pnk,nsam,ncom);
	GMM->pns = nsam;
	free(ss);

	return likelihood;
}

double** Panda_GMM_Extend_Process3(double** pnk,int nsam,int ncom,Panda_GMM *GMM)
{
	double min,*llmax,*ss,likelihood,**temppnk;
	int i,j;
	llmax=(double*)malloc(sizeof(double)*nsam);
	ss=(double*)malloc(sizeof(double)*nsam);

	temppnk = (double**)malloc(sizeof(double*)*nsam);
	for(i=0;i<nsam;i++)
		temppnk[i]=(double*)malloc(sizeof(double)*ncom);


	for(i=0;i<nsam;i++)
	{
		ss[i]=0;
		min=-DBL_MAX;
		for(j=0;j<ncom;j++)
			if (min<pnk[i][j]) min =pnk[i][j];
		llmax[i]=min;
		for(j=0;j<ncom;j++)
		{
					temppnk[i][j]=((double)exp((double)pnk[i][j]))*GMM->Coms[j].pb;
					pnk[i][j]=((double)exp(((double)(pnk[i][j]-min))))*GMM->Coms[j].pb;
					ss[i]=ss[i]+temppnk[i][j];
		}

	}

	likelihood = 0;

	//save pnk
	for(i=0;i<nsam;i++)
	{
		likelihood = likelihood + log(ss[i])+llmax[i];
		for(j=0;j<ncom;j++)
		{
			temppnk[i][j]=temppnk[i][j]/(ss[i]+DBL_MIN);
		}
	}
	free(llmax);
	free(ss);
	Panda_GMM_Util_freedouble(&pnk,nsam);
	return temppnk;
}
double *Panda_GMM_Extend_Process2(double** pnk,int nsam,int ncom,Panda_GMM *GMM)
{
	double min,*llmax,*ss,*likelihood;
	int i,j;
	llmax=(double*)malloc(sizeof(double)*nsam);
	ss=(double*)malloc(sizeof(double)*nsam);
	likelihood = (double*)malloc(sizeof(double)*nsam);
	for(i=0;i<nsam;i++)
	{
		ss[i]=0;
		min=-DBL_MAX;
		for(j=0;j<ncom;j++)
			if (min<pnk[i][j]) min =pnk[i][j];
		llmax[i]=min;
		for(j=0;j<ncom;j++)
		{
					//printf("%f   %f",pnk[i][j],exp(pnk[i][j]));
					pnk[i][j]=exp(pnk[i][j]-min)*GMM->Coms[j].pb; //range -inf -> inf
					//pnk[i][j]=exp(pnk[i][j])*GMM->Coms[j].pb; //-1 -->0
					//printf("   %f \n",pnk[i][j]);
					ss[i]=ss[i]+pnk[i][j];
		}

	}

	for(i=0;i<nsam;i++)
	{
		//likelihood[i] = log(ss[i]);//-1 -->0
		likelihood[i] = log(ss[i])+llmax[i];//range -inf -> inf
		for(j=0;j<ncom;j++)
		{
			pnk[i][j]=pnk[i][j]/ss[i];
		}
	}

	free(llmax);
	free(ss);
	return likelihood;
}
double *Panda_GMM_Extend_Process4(double** pnk,int nsam,int ncom,Panda_GMM *GMM)
{
	double min,*llmax,*ss,*likelihood;
	int i,j;
	llmax=(double*)malloc(sizeof(double)*nsam);
	ss=(double*)malloc(sizeof(double)*nsam);
	likelihood = (double*)malloc(sizeof(double)*nsam);
	for(i=0;i<nsam;i++)
	{
		ss[i]=0;
		min=-DBL_MAX;
		for(j=0;j<ncom;j++)
			if (min<pnk[i][j]) min =pnk[i][j];
		llmax[i]=min;
		for(j=0;j<ncom;j++)
		{
					//printf("%f   %f",pnk[i][j],exp(pnk[i][j]));
					pnk[i][j]=exp(pnk[i][j])*GMM->Coms[j].pb; //range -inf -> inf
					//pnk[i][j]=exp(pnk[i][j])*GMM->Coms[j].pb; //-1 -->0
					//printf("   %f \n",pnk[i][j]);
					ss[i]=ss[i]+pnk[i][j];
		}

	}

	for(i=0;i<nsam;i++)
	{
		//likelihood[i] = log(ss[i]);//-1 -->0
		likelihood[i] = 0;//range -inf -> inf
		for(j=0;j<ncom;j++)
		{
			likelihood[i]=likelihood[i]+pnk[i][j];
		}
	}

	free(llmax);
	free(ss);
	return likelihood;
}
Panda_GMM *Panda_GMM_Create(Panda_GMM_Data *data,int initK,int finalK)
{
	Panda_GMM *G,*Gmin,*Gn;
	int k;
	double resa;

    G = Panda_GMM_Init(data,initK);
    Panda_GMM_EMIterate(data,G);
	k=initK-1;
	Gmin = Panda_GMM_Util_DuplicateGMM(G);
	resa = G->rissanen;

	while (k>=finalK)
	{
        Panda_GMM_MDLReduceOrder(G);
        //qDebug()<<"reduced"<<k;
        Panda_GMM_EMIterate(data,G);
        //Panda_GMM_FreeGMM(G);
		if (resa>G->rissanen)
		{
            Panda_GMM_FreeGMM(Gmin);
			Gmin = Panda_GMM_Util_DuplicateGMM(G);
			resa = Gmin->rissanen;
		}
		k--;
	}
	Panda_GMM_FreeGMM(G);
	return Gmin;
}
double Panda_GMM_likelihood(Panda_GMM_Data *data,Panda_GMM *GMM)
{
    double *ss,*minus,*in,**pnk,likelihood;
    int K,k,i;
    ss=(double*)malloc(sizeof(double)*data->nSamples);
    memset(ss,0,sizeof(double)*data->nSamples);
    likelihood =0;

    K = GMM->NumOfMix;
    pnk = (double**)malloc(sizeof(double)*data->nSamples);
    for(int i=0;i<data->nSamples;i++)
    {
        pnk[i] = (double*)malloc(sizeof(double)*data->dim);
        memset(pnk[i],0,sizeof(double)*data->dim);
    }

    for(i=0;i<data->nSamples;i++)
    {
        for(k=0;k<K;k=k+1)
        {
            minus = Panda_GMM_Util_Minus_vecAvec(data->data[i],GMM->Coms[k].Mu,GMM->dim);
            in=Panda_GMM_Util_Multi_vecAma(minus,GMM->Coms[k].invR,GMM->dim,GMM->dim,-0.5);
            pnk[i][k]=Panda_GMM_Util_Multi_vecAvec(minus,in,GMM->dim)+GMM->Coms[k].Const;

            pnk[i][k]=((double)exp(((double)(pnk[i][k]))))*GMM->Coms[k].pb;
            ss[i]=ss[i]+pnk[i][k];

            free(minus);
            free(in);
        }
        likelihood=likelihood+log(ss[i]);
        for(k=0;k<K;k=k+1)
        {
            pnk[i][k]=pnk[i][k]/ss[i];
        }
    }

    free(ss);

    for(int i=0;i<data->nSamples;i++)
    {
        free(pnk[i]);
    }
    free(pnk);
    return likelihood/(double)data->nSamples;
}
double *Panda_GMM_Gaupdf(Panda_GMM_Data *data,Panda_GMM *GMM)
{

	double **pnk;
	double *minus,*in,*re;
	int K,k,i;

	K=GMM->NumOfMix;
	pnk = Panda_GMM_Util_Zeros(data->nSamples,GMM->NumOfMix);
	for(k=0;k<K;k=k+1)
	{
		for(i=0;i<data->nSamples;i++)
		{
			minus = Panda_GMM_Util_Minus_vecAvec(data->data[i],GMM->Coms[k].Mu,GMM->dim);
			in=Panda_GMM_Util_Multi_vecAma(minus,GMM->Coms[k].invR,GMM->dim,GMM->dim,-0.5);
			pnk[i][k]=Panda_GMM_Util_Multi_vecAvec(minus,in,GMM->dim)+GMM->Coms[k].Const;
			free(minus);
			free(in);
		}
	}

	re = Panda_GMM_Extend_Process4(pnk,data->nSamples,GMM->NumOfMix,GMM);
	Panda_GMM_Util_freedouble(&pnk,data->nSamples);
	return re;
}
/*Traingning */
void Panda_GMM_EMIterate(Panda_GMM_Data  *data,Panda_GMM *GMM)
{

    double Lc,epsilon,llnew,llold;


	Lc = 1 + GMM->dim + 0.5f * GMM->dim * (GMM->dim+1);
	epsilon = 0.01*Lc*log((double)data->nSamples*data->dim);

    llnew=Panda_GMM_EStep(data,GMM);

	while (1)
	{
		llold=llnew;
        Panda_GMM_MStep(data,GMM);
        //Panda_GMM_Util_freedouble(&pnk,data->nSamples);
        llnew=Panda_GMM_EStep(data,GMM);
        qDebug()<<"likeli " <<llnew-llold;
		if (llnew-llold < epsilon)
			break;
	}
	GMM->rissanen = -llnew+0.5*(GMM->NumOfMix*Lc-1)*log(GMM->dim * data->nSamples);
	GMM->loglikelihood=llnew;

}
double Panda_GMM_EStep(Panda_GMM_Data  *data,Panda_GMM *GMM)
{
    double *ss,*minus,*in,**pnk,likelihood;
    int K,k,i;
    ss=(double*)malloc(sizeof(double)*data->nSamples);
    memset(ss,0,sizeof(double)*data->nSamples);
    likelihood =0;

    K = GMM->NumOfMix;
    pnk = GMM->pnk;



    for(i=0;i<data->nSamples;i++)
	{
        for(k=0;k<K;k=k+1)
		{
			minus = Panda_GMM_Util_Minus_vecAvec(data->data[i],GMM->Coms[k].Mu,GMM->dim);
			in=Panda_GMM_Util_Multi_vecAma(minus,GMM->Coms[k].invR,GMM->dim,GMM->dim,-0.5);
			pnk[i][k]=Panda_GMM_Util_Multi_vecAvec(minus,in,GMM->dim)+GMM->Coms[k].Const;

            pnk[i][k]=((double)exp(((double)(pnk[i][k]))))*GMM->Coms[k].pb;
            ss[i]=ss[i]+pnk[i][k];

			free(minus);
			free(in);
		}
        likelihood=likelihood+log(ss[i]);
        for(k=0;k<K;k=k+1)
        {
            pnk[i][k]=pnk[i][k]/ss[i];
        }    
	}    
    GMM->pns = data->nSamples;  
    free(ss);
    return likelihood;
}
void Panda_GMM_MStep(Panda_GMM_Data  *data,Panda_GMM *GMM)
{
	int i,j,r,s,l;
	double **pnk,temp;
    pnk=GMM->pnk;


	for (i=0;i<GMM->NumOfMix;i++)
	{
		//----------
		//Update N, pb
		temp=0;
		for(j=0;j<data->nSamples;j++) temp=temp+pnk[j][i];
		GMM->Coms[i].N= temp;
		//----------
		GMM->Coms[i].pb= temp;
		//-----
		//Update Mean

		for(j=0;j<GMM->dim;j++)
		{
			temp=0;
			for(r=0;r<data->nSamples;r++)
				temp=temp+data->data[r][j]*pnk[r][i];
			temp=temp/GMM->Coms[i].N;
			GMM->Coms[i].Mu[j]=temp;
		}
		//Updata Covar matrix
		for(r=0;r<GMM->dim;r++)
			for(s=r;s<data->dim;s++)
			{
				temp=0;
				for(l=0;l<data->nSamples;l++)
				{
					temp=temp+   (data->data[l][r]-GMM->Coms[i].Mu[r])    *  ((data->data[l][s]-GMM->Coms[i].Mu[s])*pnk[l][i]);
				}
				temp=temp/GMM->Coms[i].N;
				GMM->Coms[i].R[r][s]=temp;
				GMM->Coms[i].R[s][r]=temp;
				if (r==s) GMM->Coms[i].R[r][s]=GMM->Coms[i].R[r][s]+GMM->Rmin;
			}

	}
	Panda_GMM_Normalize(GMM);	
}
Panda_GMM_Component *Panda_GMM_Component_Combine(Panda_GMM_Component *Com1,Panda_GMM_Component *Com2,int dim)
{
	Panda_GMM_Component *re;
	double wt1,wt2;
	int i,j;
	re=Panda_GMM_Util_DuplicateGMMComponent(Com1,dim);

	wt1 = Com1->N / (Com1->N+Com2->N);
	wt2 = 1.0f-wt1;
	for(i=0;i<dim;i++)
		re->Mu[i]=Com1->Mu[i]*wt1+Com2->Mu[i]*wt2;

	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			{
				re->R[i][j]= wt1*(Com1->R[i][j] + (re->Mu[i]-Com1->Mu[i])*(re->Mu[j]-Com1->Mu[j]))
						+ wt2*(Com2->R[i][j] + (re->Mu[i]-Com2->Mu[i])*(re->Mu[j]-Com2->Mu[j]));
			}
	for(i=0;i<dim;i++)
			for(j=0;j<dim;j++) { re->invR[i][j]=re->R[i][j] ;}
	invert(re->invR,dim);

	re->pb = Com1->pb + Com2->pb;
	re->N = Com1->N + Com2->N;
	re->Const = -(((double)dim)*log(2.0f*3.1416f) + ((double)log(matrix_det(re->R,dim))))/2;
	return re;
}
double Panda_GMM_Component_Distance(Panda_GMM_Component *Com1,Panda_GMM_Component *Com2,int dim)
{
	double re;
	Panda_GMM_Component *Com3;
	Com3=Panda_GMM_Component_Combine(Com1,Com2,dim);
	re= Com1->Const*Com1->N + Com2->Const*Com2->N - Com3->Const*Com3->N;
	Panda_GMM_FreeComponent(Com3,dim);
	return re;
}
Panda_GMM *Panda_GMM_MDLReduceOrder(Panda_GMM *GMM)
{
	int i,j,mink1,mink2;
	double dmin,temp;
	Panda_GMM_Component *com3;
	Panda_GMM *re;
	dmin=0;
	for(i=0;i<GMM->NumOfMix;i++)
		for(j=i+1;j<GMM->NumOfMix;j++)
		{
			temp=Panda_GMM_Component_Distance(&(GMM->Coms[i]),&(GMM->Coms[j]),GMM->dim);
			if ((i==0 && j==1)||dmin>temp)
			{
				dmin=temp;
				mink1=i;mink2=j;
			}
		}
	///////////////////
//	re = (Panda_GMM*)malloc(sizeof(Panda_GMM));
//	re->NumOfMix = GMM->NumOfMix-1;
//	re->dim=GMM->dim;
//	re->Coms = (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component)*re->NumOfMix);
//	for(i=0;i<re->NumOfMix;i++)
//		Panda_GMM_Component_Init(&(re->Coms[i]),re->dim);
//	re->loglikelihood=GMM->loglikelihood;
//	re->Rmin=GMM->Rmin;
//	re->rissanen=GMM->rissanen;
//    re->pnk = (double**) malloc(sizeof(double*)*GMM->pns);
//    re->pns = GMM->pns;
//    for(i=0;i<GMM->pns;i++)
//    {
//        re->pnk[i] = (double*) malloc(sizeof(double)*re->NumOfMix);
//    }

    GMM->NumOfMix = GMM->NumOfMix-1;
    com3 = Panda_GMM_Component_Combine(GMM->Coms+mink1,GMM->Coms + mink2,GMM->dim);
//	for(i=0;i<mink1;i++)
//		Panda_GMM_Component_Copy(&(re->Coms[i]),&(GMM->Coms[i]),GMM->dim);
//	Panda_GMM_Component_Copy(&(re->Coms[i]),com3,GMM->dim);
//	for(i=mink1+1;i<mink2;i++)
//		Panda_GMM_Component_Copy(&(re->Coms[i]),&(GMM->Coms[i]),GMM->dim);
//	for(i=mink2+1;i<re->NumOfMix+1;i++)
//			Panda_GMM_Component_Copy(&(re->Coms[i-1]),&(GMM->Coms[i]),GMM->dim);
    Panda_GMM_FreeComponent(GMM->Coms + mink2,GMM->dim);
    Panda_GMM_FreeComponent(GMM->Coms + mink1,GMM->dim);

    memcpy(GMM->Coms+mink1,com3,sizeof(Panda_GMM_Component));

    for(i=mink2+1;i<GMM->NumOfMix+1;i++)
        memcpy(GMM->Coms+i-1,GMM->Coms+i,sizeof(Panda_GMM_Component));

    //Panda_GMM_FreeComponent(GMM->Coms + GMM->NumOfMix-1,GMM->dim);
    free(com3);
    Panda_GMM_Normalize(GMM);



    return GMM;
}
/*Increamental Learning*/
Panda_GMM *Panda_GMM_IncreamentalLearning(Panda_GMM *GMM,Panda_GMM_Data *data)
{
	double loglik_threshod = 1e-10,loglik,loglik_old = - DBL_MAX, *E0,*E,**pnk,su,**covtmp,*li;
	int i,j,z,h;
	Panda_GMM *re;
	re = Panda_GMM_Util_DuplicateGMM(GMM);

	E0 = Panda_GMM_IncreamentalLearning_SumPnk(GMM->pnk,GMM->pns,GMM->NumOfMix);
	pnk = NULL;
	E = NULL;

	int nbStep =0;
	while (1)
	{
		if (pnk != NULL)
		{
			Panda_GMM_Util_freedouble(&pnk,data->nSamples);
		}
		pnk = Panda_GMM_IncreamentalLearning_pnk(re,data);

		if (E!=NULL) free(E);
		E = Panda_GMM_IncreamentalLearning_SumPnk(pnk,data->nSamples,GMM->NumOfMix);

		for(i=0;i<GMM->NumOfMix;i++)
		{
			//update prors
			re->Coms[i].pb = (E0[i]+E[i]) / ((double) (data->nSamples + GMM->pns));
			//update mean

			for(j=0;j<GMM->dim;j++)
			{
				su=0;
				for(z=0;z<data->nSamples;z++)
					su=su+data->data[z][j]*pnk[z][i];
				re->Coms[i].Mu[j]=(GMM->Coms[i].Mu[j]*E0[i]+su)/(E0[i]+E[i]);
			}
			//Update Covariance
			covtmp = Panda_GMM_Util_Zeros(GMM->dim,GMM->dim);
			for(j=0;j<data->nSamples;j++)
			{
				for(z=0;z<GMM->dim;z++)
					for(h=0;h<GMM->dim;h++)
						covtmp[z][h]=covtmp[z][h]+(data->data[j][z]-re->Coms[i].Mu[z])*(data->data[j][h]-re->Coms[i].Mu[h])*pnk[j][i];
			}
			for(z=0;z<GMM->dim;z++)
				for(h=0;h<GMM->dim;h++)
				{
					re->Coms[i].R[z][h]=((GMM->Coms[i].R[z][h]+
							(GMM->Coms[i].Mu[z]-re->Coms[i].Mu[z])*(GMM->Coms[i].Mu[h]-re->Coms[i].Mu[h]))*E0[i]+covtmp[z][h])/(E0[i]+E[i]);
				}
			Panda_GMM_Util_freedouble(&covtmp,GMM->dim);
		}
		Panda_GMM_Normalize_withoutPb(re);
		loglik=0;
        //li=Panda_GMM_likelihood(data,re);
		for(j=0;j<data->nSamples;j++)
			loglik = loglik + li[j];
		free(li);
		loglik = loglik / (double)data->nSamples;

		if (abs(loglik / loglik_old -1)<loglik_threshod) break;

		loglik_old =loglik;

		nbStep ++;


	}


	Panda_GMM_Normalize(re);
	free(E0);
	Panda_GMM_IncreamentalLearning_PlusPnk(re,pnk,data->nSamples);
	if (pnk != NULL)
	{
		Panda_GMM_Util_freedouble(&pnk,data->nSamples);
	}
	return re;
}
double** Panda_GMM_IncreamentalLearning_pnk(Panda_GMM *GMM,Panda_GMM_Data *data)
{
	double *minus,*in,**pnk;

	int K,k,i;

	K=GMM->NumOfMix;
	pnk = Panda_GMM_Util_Zeros(data->nSamples,GMM->NumOfMix);

	for(k=0;k<K;k=k+1)
	{
		for(i=0;i<data->nSamples;i++)
		{
			minus = Panda_GMM_Util_Minus_vecAvec(data->data[i],GMM->Coms[k].Mu,GMM->dim);
			in=Panda_GMM_Util_Multi_vecAma(minus,GMM->Coms[k].invR,GMM->dim,GMM->dim,-0.5);
			pnk[i][k]=Panda_GMM_Util_Multi_vecAvec(minus,in,GMM->dim)+GMM->Coms[k].Const;
			free(minus);
			free(in);
		}
	}
	pnk=Panda_GMM_Extend_Process3(pnk,data->nSamples,GMM->NumOfMix,GMM);
	return pnk;
}
void Panda_GMM_IncreamentalLearning_PlusPnk(Panda_GMM *GMM,double **pnk,int nsam)
{
	int tnsam,i,j,z;
	double **re;
	tnsam = nsam + GMM->pns;

	re=(double**)malloc(sizeof(double*)*tnsam);
	for(i=0;i<tnsam;i++) re[i]=(double*)malloc(sizeof(double)*GMM->NumOfMix);

	for(i=0;i<GMM->pns;i++)
	{
		for(j=0;j<GMM->NumOfMix;j++)
			re[i][j] = GMM->pnk[i][j];
	}
	for(z=0;z<nsam;z++)
	{
		for(j=0;j<GMM->NumOfMix;j++)
		{
			re[i][j]=pnk[z][j];
		}
		i++;
	}
	Panda_GMM_Util_freedouble(&(GMM->pnk),GMM->pns);
	GMM->pnk = re;
	GMM->pns = tnsam;
}
double* Panda_GMM_IncreamentalLearning_SumPnk(double **pnk,int nsam,int NumMix)
{
	double *re,su;
	int i,j;
	re = (double*)malloc(sizeof(double)*NumMix);
	for(i=0;i<NumMix;i++)
	{
		su=0;
		for(j=0;j<nsam;j++)
		{
			su =su + pnk[j][i];
		}
		re[i]=su;
	}
	return re;
}
//save file
void Panda_GMM_Save(char *filepath,Panda_GMM* GMM)
{
	int i,j;
	FILE *f;
	f = fopen(filepath,"wb");
	fwrite(&(GMM->NumOfMix),sizeof(int),1,f);
	fwrite(&(GMM->dim),sizeof(int),1,f);
	fwrite(&(GMM->Rmin),sizeof(double),1,f);
	fwrite(&(GMM->rissanen),sizeof(double),1,f);
	fwrite(&(GMM->loglikelihood),sizeof(double),1,f);
	fwrite(&(GMM->pns),sizeof(int),1,f);

	for(i=0 ;i<GMM->NumOfMix;i++)
	{
		fwrite(&(GMM->Coms[i].N),sizeof(double),1,f);
		fwrite(&(GMM->Coms[i].pb),sizeof(double),1,f);
		fwrite(&(GMM->Coms[i].Const),sizeof(double),1,f);
		fwrite((GMM->Coms[i].Mu),sizeof(double),GMM->dim,f);

		for(j=0;j<GMM->dim;j++)
			fwrite((GMM->Coms[i].R[j]),sizeof(double),GMM->dim,f);
		for(j=0;j<GMM->dim;j++)
					fwrite((GMM->Coms[i].invR[j]),sizeof(double),GMM->dim,f);
	}

	if (GMM->pns!=0)
	{
		for(i=0;i<GMM->pns;i++)
			fwrite((GMM->pnk[i]),sizeof(double),GMM->NumOfMix,f);
	}    
    fwrite((GMM->name),sizeof(char),100,f);
	fclose(f);
}
Panda_GMM* Panda_GMM_Read(char *filepath)
{
	Panda_GMM* GMM;
	int i,j;
	FILE *f;
	f = fopen(filepath,"rb");
	GMM = (Panda_GMM*)malloc(sizeof(Panda_GMM));
	fread(&(GMM->NumOfMix),sizeof(int),1,f);
	fread(&(GMM->dim),sizeof(int),1,f);
	fread(&(GMM->Rmin),sizeof(double),1,f);
	fread(&(GMM->rissanen),sizeof(double),1,f);
	fread(&(GMM->loglikelihood),sizeof(double),1,f);
	fread(&(GMM->pns),sizeof(int),1,f);

	GMM->Coms = (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component)*GMM->NumOfMix);
	for(i=0 ;i<GMM->NumOfMix;i++)
	{
		fread(&(GMM->Coms[i].N),sizeof(double),1,f);
		fread(&(GMM->Coms[i].pb),sizeof(double),1,f);
		fread(&(GMM->Coms[i].Const),sizeof(double),1,f);

		GMM->Coms[i].Mu=(double*)malloc(sizeof(double)*GMM->dim);
		fread((GMM->Coms[i].Mu),sizeof(double),GMM->dim,f);

		GMM->Coms[i].R = (double**)malloc(sizeof(double*)*GMM->dim);
		GMM->Coms[i].invR = (double**)malloc(sizeof(double*)*GMM->dim);

		for(j=0;j<GMM->dim;j++)
		{
			GMM->Coms[i].R[j] = (double*)malloc(sizeof(double)*GMM->dim);
			fread((GMM->Coms[i].R[j]),sizeof(double),GMM->dim,f);
		}
		for(j=0;j<GMM->dim;j++)
		{
			GMM->Coms[i].invR[j] = (double*)malloc(sizeof(double)*GMM->dim);
			fread((GMM->Coms[i].invR[j]),sizeof(double),GMM->dim,f);
		}
	}

	if (GMM->pns!=0)
	{
		GMM->pnk = (double**) malloc(sizeof(double*)*GMM->pns);
		for(i=0;i<GMM->pns;i++)
		{
			GMM->pnk[i]=(double*)malloc(sizeof(double)*GMM->NumOfMix);
			fread((GMM->pnk[i]),sizeof(double),GMM->NumOfMix,f);
		}
	}

    fread((GMM->name),sizeof(char),100,f);
	fclose(f);
	return GMM;
}
Panda_GMM* Panda_GMM_Load_backgroundModel(char *filepath)
{
    Panda_GMM* GMM;
    int i,j;
    FILE *f;
    f = fopen(filepath,"rb");
    GMM = (Panda_GMM*)malloc(sizeof(Panda_GMM));
    fread(&(GMM->NumOfMix),sizeof(int),1,f);
    fread(&(GMM->dim),sizeof(int),1,f);
    fread(&(GMM->Rmin),sizeof(double),1,f);
    fread(&(GMM->rissanen),sizeof(double),1,f);
    fread(&(GMM->loglikelihood),sizeof(double),1,f);
    fread(&(GMM->pns),sizeof(int),1,f);

    GMM->Coms = (Panda_GMM_Component*)malloc(sizeof(Panda_GMM_Component)*GMM->NumOfMix);
    for(i=0 ;i<GMM->NumOfMix;i++)
    {
        fread(&(GMM->Coms[i].N),sizeof(double),1,f);
        fread(&(GMM->Coms[i].pb),sizeof(double),1,f);
        fread(&(GMM->Coms[i].Const),sizeof(double),1,f);

        GMM->Coms[i].Mu=(double*)malloc(sizeof(double)*GMM->dim);
        fread((GMM->Coms[i].Mu),sizeof(double),GMM->dim,f);

        GMM->Coms[i].R = (double**)malloc(sizeof(double*)*GMM->dim);
        GMM->Coms[i].invR = (double**)malloc(sizeof(double*)*GMM->dim);

        for(j=0;j<GMM->dim;j++)
        {
            GMM->Coms[i].R[j] = (double*)malloc(sizeof(double)*GMM->dim);
            fread((GMM->Coms[i].R[j]),sizeof(double),GMM->dim,f);
        }
        for(j=0;j<GMM->dim;j++)
        {
            GMM->Coms[i].invR[j] = (double*)malloc(sizeof(double)*GMM->dim);
            fread((GMM->Coms[i].invR[j]),sizeof(double),GMM->dim,f);
        }
    }

    if (GMM->pns!=0)
    {
        GMM->pnk = (double**) malloc(sizeof(double*)*GMM->pns);
        for(i=0;i<GMM->pns;i++)
        {
            GMM->pnk[i]=(double*)malloc(sizeof(double)*GMM->NumOfMix);
            fread((GMM->pnk[i]),sizeof(double),GMM->NumOfMix,f);
        }
    }
    fclose(f);
    return GMM;
}
