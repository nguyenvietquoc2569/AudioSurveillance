
#include "GMManalyser.h"
#include "utils.h"
#include <qmath.h>
#include <qmetatype.h>
#include <QAudioFormat>
#include <QThread>
#include "mfccLib.h"
#include "GMManalyser.h"
#include <QDir>

const   int initmaxK = 7;
const   int initminK = 1;

GMMAnalyserThread::GMMAnalyserThread(QObject *parent,Panda_mfccParameter *para,QString _Background_folder)
    :   QObject(parent)
    ,   m_thread(new QThread(this))
    ,   background_folder(_Background_folder)

{
    // moveToThread() cannot be called on a QObject with a parent
    mf_format = para;
    setParent(0);
    m_thread->setPriority(QThread::HighPriority);
    moveToThread(m_thread);
    m_thread->start();


}

GMMAnalyserThread::~GMMAnalyserThread()
{
}

void GMMAnalyserThread::calculateGMM_background(const PBufferFeature &buffer,qint64 pos,qint64 length)
{
    qDebug() << "GMMAnalyserThread::calculate Init GMM" << QThread::currentThread();

    Panda_GMM *GMM;
    Panda_GMM_Data *GMM_Data;
    qreal result;


    //------------------
    //------------------
    QList<double*> temp(buffer.list);
    QList<double*>::iterator i;
    double** re = (double**)malloc(buffer.list.count()*sizeof(double*));
    int dem=0;
    for(i=temp.begin();i!=temp.end();++i)
    {
        re[dem]=*i;
        dem++;
    }
    //GMM_Data = new Panda_GMM_Data(re,mf_format->M,dem);
    GMM_Data = (Panda_GMM_Data*)malloc(sizeof(Panda_GMM_Data));
    GMM_Data->data=re;
    GMM_Data->dim=mf_format->M;
    GMM_Data->nSamples=dem;
    //--------------------
    //--------------------

    QDir dir(background_folder);
    dir.setNameFilters(QStringList()<<"*.GMM");
    dir.setFilter(QDir::Files);

    Panda_GMM *GMMtemp;
    for(int j=0;j<dir.entryList().count();j++)
    {
        if (j==0)
        {
            QString qpath = dir.absoluteFilePath(dir.entryList()[j]);
            QByteArray Apath = qpath.toLocal8Bit();
            GMM = Panda_GMM_Load_backgroundModel(Apath.data());

            GMM->likelihood = Panda_GMM_likelihood(GMM_Data,GMM);


        }
        else
        {
            QString qpath = dir.absoluteFilePath(dir.entryList()[j]);
            QByteArray Apath = qpath.toLocal8Bit();
            GMMtemp = Panda_GMM_Load_backgroundModel(Apath.data());
            GMMtemp->likelihood = Panda_GMM_likelihood(GMM_Data,GMM);

            if (GMMtemp->likelihood > GMM->likelihood)
            {
                qDebug()<<"the Chossen model :"<<dir.absoluteFilePath(dir.entryList()[j]);
                free(GMM);
                GMM=GMMtemp;
            }

        }

    }

    qDebug()<<"GMMAnalyserThread::calculateGMM :  finish training with  "<<dem <<GMM_Data->nSamples <<" samples";
    qDebug()<<"Gmm : number component " <<GMM->NumOfMix;
    qDebug()<<"initial likelihood " <<result;

    free(GMM_Data->data);
    free(GMM_Data);

    emit calculationComplete_background(GMM,pos,length);
}

//void GMMAnalyserThread::calculateGMM_background(const PBufferFeature &buffer,qint64 pos,qint64 length)
//{
//    qDebug() << "GMMAnalyserThread::calculate Init GMM" << QThread::currentThread();

//    Panda_GMM *GMM;
//    Panda_GMM_Data *GMM_Data;
//    qreal result;


//    //------------------
//    //------------------
//    QList<double*> temp(buffer.list);
//    QList<double*>::iterator i;
//    double** re = (double**)malloc(buffer.list.count()*sizeof(double*));
//    int dem=0;
//    for(i=temp.begin();i!=temp.end();++i)
//    {
//        re[dem]=*i;
//        dem++;
//    }
//    GMM_Data = new Panda_GMM_Data(re,mf_format->M,dem);
//    //--------------------
//    //--------------------

//    QDir dir(Background_Folder);
//    dir.setNameFilters(QStringList()<<"*.GMM");
//    dir.setFilter(QDir::Files);;

//    Panda_GMM *GMMtemp;

//    GMM = Panda_GMM_Create(GMM_Data,4,1);
//    double* llar = Panda_GMM_likelihood(GMM_Data,GMM);
//    GMM->likelihood = MeanLL(llar,GMM_Data->nSamples);
//    Panda_GMM_Save("d:\\Quoc1.GMM",GMM);

//    free(llar);

//    qDebug()<<"GMMAnalyserThread::calculateGMM :  finish training with  "<<dem <<GMM_Data->nSamples <<" samples";
//    qDebug()<<"Gmm : number component " <<GMM->NumOfMix;
//    qDebug()<<"initial likelihood " <<result;

//    free(GMM_Data->data);
//    free(GMM_Data);

//    emit calculationComplete_background(GMM,pos,length);
//}

void GMMAnalyserThread::calculateGMM(const PBufferFeature &buffer,qint64 pos,qint64 length)
{
    qDebug() << "GMMAnalyserThread::calculate Init GMM"
                           << QThread::currentThread();

    Panda_GMM *GMM;
    Panda_GMM_Data *GMM_Data;

    QList<double*> temp(buffer.list);
    QList<double*>::iterator i;
    double** re = (double**)malloc(buffer.list.count()*sizeof(double*));
    int dem=0;
    for(i=temp.begin();i!=temp.end();++i)
    {
        re[dem]=*i;
        dem++;
    }


    GMM_Data = (Panda_GMM_Data*)malloc(sizeof(Panda_GMM_Data));
    GMM_Data->data=re;
    GMM_Data->dim=mf_format->M;
    GMM_Data->nSamples=dem;



    //GMM = Panda_GMM_Create(GMM_Data,initmaxK,initminK);
    GMM = Panda_GMM_Read("D:\\GMMBackground.GMM");
    qDebug()<<"GMMAnalyserThread::calculateGMM :  finish training with  "<<dem <<GMM_Data->nSamples <<" samples";
    qDebug()<<"Gmm : number component " <<GMM->NumOfMix;
    qreal result = Panda_GMM_likelihood(GMM_Data,GMM);;
    qDebug()<<"initial likelihood " <<result;
    GMM->likelihood = result;
    free(GMM_Data->data);
    free(GMM_Data);

    emit calculationComplete_background(GMM,pos,length);
}

//=============================================================================
// GMMAnalyser
//=============================================================================

GMMAnalyser::GMMAnalyser(QObject *parent,Panda_mfccParameter *para,QString folderlink)
    :   QObject(parent)
    ,   m_thread(new GMMAnalyserThread(this,para,folderlink))
    ,   m_state(Idle)
#ifdef DUMP_GMMANALYSER
    ,   m_count(0)
#endif
{
    Background_folder = folderlink;
    CHECKED_CONNECT(m_thread, SIGNAL(calculationComplete(Panda_GMM*,qint64,qint64)),
                    this, SLOT(calculationComplete(Panda_GMM*,qint64,qint64)));
    CHECKED_CONNECT(m_thread, SIGNAL(calculationComplete_background(Panda_GMM*,qint64,qint64)),
                    this, SLOT(calculationComplete_background(Panda_GMM*,qint64,qint64)));
}

GMMAnalyser::~GMMAnalyser()
{

}

//-----------------------------------------------------------------------------
// Public functions
//-----------------------------------------------------------------------------


void GMMAnalyser::calculate(const PBufferFeature &buffer,
                         qint64 pos,
                         qint64 length)
{

    // QThread::currentThread is marked 'for internal use only', but
    // we're only using it for debug output here, so it's probably OK :)
    qDebug() << "GMMAnalyser::calculate"
                           << QThread::currentThread();


    //if (isReady()) {
        m_state = Busy;

        // Invoke GMMAnalyserThread::calculateGMM using QMetaObject.  If
        // m_thread is in a different thread from the current thread, the
        // calculation will be done in the child thread.
        // Once the calculation is finished, a calculationChanged signal will be
        // emitted by m_thread.
        const bool b = QMetaObject::invokeMethod(m_thread, "calculateGMM",
                                  Qt::AutoConnection,
                                  Q_ARG(PBufferFeature, buffer),
                                  Q_ARG(qint64, pos),
                                  Q_ARG(qint64, length));
        Q_ASSERT(b);
        Q_UNUSED(b) // suppress warnings in release builds
    //}
}
void GMMAnalyser::calculate_backgournd(const PBufferFeature &buffer,
                         qint64 pos,
                         qint64 length)
{

    // QThread::currentThread is marked 'for internal use only', but
    // we're only using it for debug output here, so it's probably OK :)
    qDebug() << "GMMAnalyser::calculate"
                           << QThread::currentThread();


    //if (isReady()) {
        m_state = Busy;

        // Invoke GMMAnalyserThread::calculateGMM using QMetaObject.  If
        // m_thread is in a different thread from the current thread, the
        // calculation will be done in the child thread.
        // Once the calculation is finished, a calculationChanged signal will be
        // emitted by m_thread.
        const bool b = QMetaObject::invokeMethod(m_thread, "calculateGMM_background",
                                  Qt::AutoConnection,
                                  Q_ARG(PBufferFeature, buffer),
                                  Q_ARG(qint64, pos),
                                  Q_ARG(qint64, length));
        Q_ASSERT(b);
        Q_UNUSED(b) // suppress warnings in release builds
    //}
}

bool GMMAnalyser::isReady() const
{
    return (Idle == m_state);
}

void GMMAnalyser::cancelCalculation()
{
    if (Busy == m_state)
        m_state = Cancelled;
}


//-----------------------------------------------------------------------------
// Private slots
//-----------------------------------------------------------------------------

void GMMAnalyser::calculationComplete(Panda_GMM *GMM,qint64 pos,qint64 length)
{
    //Q_ASSERT(Idle != m_state);
    //if (Busy == m_state)
    this->m_GMM = GMM;
    emit GMMInitComplete(GMM,pos,length);
    //m_state = Idle;
}
void GMMAnalyser::calculationComplete_background(Panda_GMM *GMM,qint64 pos,qint64 length)
{
    //Q_ASSERT(Idle != m_state);
    //if (Busy == m_state)
    this->m_GMM = GMM;
    emit GMMInitComplete_background(GMM,pos,length);
    //m_state = Idle;
}
