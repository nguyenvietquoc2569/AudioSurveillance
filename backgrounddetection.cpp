#include "backgrounddetection.h"
#include "utils.h"
#include <QDebug>
#include <QDir>
#include <float.h>


BackgroundDetection::BackgroundDetection(QObject *parent) :
    QObject(parent)
{
    mf_format = 0;
    GMM = 0;
    GMMCount = 0;

    this->threshold=50;
}

void BackgroundDetection::ReadClassifierFromFolder(QString folder)
{
    QDir dir(folder);
    dir.setNameFilters(QStringList()<<"*.scl");
    dir.setFilter(QDir::Files);
    GMMSound = (Panda_GMM**)malloc(sizeof(Panda_GMM*)*dir.entryList().count());
    qDebug()<<"Reading classifier";
    for(int j=0;j<dir.entryList().count();j++)
    {

            QString qpath = dir.absoluteFilePath(dir.entryList()[j]);
            QByteArray Apath = qpath.toLocal8Bit();
            qDebug()<<qpath;
            GMMSound[GMMCount] = Panda_GMM_Read(Apath.data());
            GMMCount++;
    }
}

void BackgroundDetection::GMMChange(Panda_GMM *GMM,qint64 pos,qint64 length)
{
    this->GMM=GMM;
}

void BackgroundDetection::FeaturebufferChange(double **arr,int num,qint64 pos,qint64 length)
{
    //Panda_GMM_Data  *data = new Panda_GMM_Data(arr,mf_format->M,num);
    //double* llar = Panda_GMM_likelihood(data,GMM);

    Panda_GMM_Data  *data;
    data = (Panda_GMM_Data*)malloc(sizeof(Panda_GMM_Data));
    data->data=arr;
    data->dim=mf_format->M;
    data->nSamples=num;

    qreal result = Panda_GMM_likelihood(data,GMM);
    qDebug()<<"log likelihood :" <<result << "distance" << qAbs(result-GMM->likelihood)<<"number of data "<<num;



    if (qAbs(result-GMM->likelihood)>threshold)
    {
        QString classname="";
        double likeli = -DBL_MAX;
        for (int i = 0; i < GMMCount; ++i)
        {
            double tam = Panda_GMM_likelihood(data,GMMSound[i]);
            if (tam>likeli)
            {
                likeli=tam;
                classname = QString(GMMSound[i]->name);
            }
        }
        emit ForegroundDetected(pos,length,classname);
    }
    else
    {
        emit BackgroundDetected(pos,length);
    }
    free(data);free(arr);
}
