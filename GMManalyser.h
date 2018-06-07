#ifndef GMMANALYSER_H
#define GMMANALYSER_H
#define CHECKED_CONNECT(source, signal, receiver, slot) \
    if (!connect(source, signal, receiver, slot)) \
        qt_assert_x(Q_FUNC_INFO, "CHECKED_CONNECT failed", __FILE__, __LINE__);


#include <QByteArray>
#include <QObject>
#include <QVector>
#include "mfccLib.h"
#include "Panda_GMM.h"
#include "Panda_GMM_Util.h"
#include "pbufferfeature.h"

#ifdef DUMP_GMMANALYSER
#include <QDir>
#include <QFile>
#include <QTextStream>
#endif


#ifndef DISABLE_FFT
#endif

QT_FORWARD_DECLARE_CLASS(QAudioFormat)
QT_FORWARD_DECLARE_CLASS(QThread)

class FFTRealWrapper;

class GMMAnalyserThreadPrivate;

/**
 * Implementation of the GMM analysis which can be run in a
 * separate thread.
 */
class GMMAnalyserThread : public QObject
{
    Q_OBJECT

public:
    GMMAnalyserThread(QObject *parent,Panda_mfccParameter *para,QString _Background_folder);
    ~GMMAnalyserThread();
    Panda_mfccParameter *mf_format;    
    QString background_folder;


public slots:

    void calculateGMM(const PBufferFeature &buffer,
                           qint64 pos,
                           qint64 length);
    void calculateGMM_background(const PBufferFeature &buffer, qint64 pos, qint64 length);
signals:
    void calculationComplete(Panda_GMM *GMM,qint64 pos,qint64 length);
    void calculationComplete_background(Panda_GMM *GMM,qint64 pos,qint64 length);



private:

    QThread*                                    m_thread;

};

/**
 * Class which performs frequency GMM analysis on a window of
 * audio samples, provided to it by the Engine.
 */
class GMMAnalyser : public QObject
{
    Q_OBJECT

public:
    GMMAnalyser(QObject *parent = 0,Panda_mfccParameter *para=0,QString Folderlink="");
    ~GMMAnalyser();
    Panda_mfccParameter *mf_format;
    Panda_GMM *m_GMM;
    QList<Panda_GMM> SoundModel;

    QString Background_folder;



public:

    void calculate(const PBufferFeature &buffer, qint64 pos, qint64 lenth);
    bool isReady() const;
    void cancelCalculation();    
    void calculate_backgournd(const PBufferFeature &buffer, qint64 pos, qint64 length);
signals:
    void GMMChanged(Panda_GMM *GMM,qint64,qint64);
    void GMMInitComplete(Panda_GMM *GMM,qint64,qint64);
    void GMMInitComplete_background(Panda_GMM *GMM,qint64,qint64);

private slots:
    void calculationComplete(Panda_GMM *GMM,qint64,qint64);
    void calculationComplete_background(Panda_GMM *GMM,qint64,qint64);

private:
    void calculateWindow();

private:
    GMMAnalyserThread*    m_thread;

    enum State {
        Idle,
        Busy,
        Cancelled
    };

    State              m_state;
};

#endif // GMMANALYSER_H

