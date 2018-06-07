#ifndef BACKGROUNDDETECTION_H
#define BACKGROUNDDETECTION_H

#include <QObject>
#include "Panda_GMM.h"
#include "mfccLib.h"
#include "pbufferfeature.h"

class BackgroundDetection : public QObject
{
    Q_OBJECT
public:
    explicit BackgroundDetection(QObject *parent = 0);
    qreal threshold;
    Panda_GMM *GMM,**GMMSound;
    Panda_mfccParameter *mf_format;
    int GMMCount;

    void ReadClassifierFromFolder(QString folder);
    
signals:
    void ForegroundDetected(qint64 pos,qint64 length,QString classname);
    void BackgroundDetected(qint64 pos,qint64 length);
    
public slots:
    void GMMChange(Panda_GMM *GMM, qint64 pos,qint64 length);
    void FeaturebufferChange(double **arr,int num,qint64 pos,qint64 length);

private:



};

#endif // BACKGROUNDDETECTION_H
