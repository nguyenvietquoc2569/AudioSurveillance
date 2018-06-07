

#ifndef SSLANALYSER_H
#define SSLANALYSER_H
#define CHECKED_CONNECT(source, signal, receiver, slot) \
    if (!connect(source, signal, receiver, slot)) \
        qt_assert_x(Q_FUNC_INFO, "CHECKED_CONNECT failed", __FILE__, __LINE__);


#include <QByteArray>
#include <QObject>
#include <QVector>
#include "mfccLib.h"
#include "fftssl.h"

#ifdef DUMP_SSLANALYSER
#include <QDir>
#include <QFile>
#include <QTextStream>
#endif


#ifndef DISABLE_FFT
#endif

QT_FORWARD_DECLARE_CLASS(QAudioFormat)
QT_FORWARD_DECLARE_CLASS(QThread)

class FFTRealWrapper;

class SSLAnalyserThreadPrivate;

/**
 * Implementation of the SSL analysis which can be run in a
 * separate thread.
 */
class SSLAnalyserThread : public QObject
{
    Q_OBJECT

public:
    SSLAnalyserThread(QObject *parent, double baseline,double soundspeed, int nDelays, double beta, qint64 FFTSize,qint64 shiftFrame,qint64 SampleRate);
    ~SSLAnalyserThread();
    double CalculateForeachPart(const QByteArray &buffer1,const QByteArray &buffer2,qint64 pos,qint64 length);
    void calculatemCoinc(double *leftFFT, double *rightFFT);
    void initAz();
    FFTSSL *mFFt;//
    qint64 mNDelays;//
    qint64 mFFTSize;//
    qint64 mShiftFrame;
    qint64 mFFTSize_2;//
    double* mDelayFactors;//
    double* mCoinMap;//
    double mBeta;//
    double* mWindow;//
    double* mLeftIn,*mRightIn,* mLeftIn1,*mRightIn1;//
    double *mAzVals;
    double msoundspeed;
    qint64 mSamplerate;

public slots:

    void calculateSSL(const QByteArray &buffer,
                      const QByteArray &buffer1,
                      const QByteArray &buffer2,
                      const QByteArray &buffer3,
                   qint64 pos,
                   qint64 length);

signals:
    void calculationComplete(double angel1,double angel2,qint64,qint64);



private:

    QThread*                                    m_thread;

};

/**
 * Class which performs frequency SSL analysis on a window of
 * audio samples, provided to it by the Engine.
 */
class SSLAnalyser : public QObject
{
    Q_OBJECT
public:
    SSLAnalyser(QObject *parent, double baseline,double soundspeed, int nDelays, double beta, qint64 FFTSize,qint64 shiftFrame,qint64 SampleRate);
    ~SSLAnalyser();
    double* mLeftIn,*mRightIn,* mLeftIn1,*mRightIn1;







public:
    void calculate(const QByteArray &buffer,
                   const QByteArray &buffer1,
                   const QByteArray &buffer2,
                   const QByteArray &buffer3,
                qint64 pos,
                qint64 length);
    bool isReady() const;

    void cancelCalculation();
signals:
    void SSLChanged(double angel1,double angel2,qint64,qint64);
private slots:
    void calculationComplete(double angel1,double angel2,qint64,qint64);
private:
    void calculateWindow();
private:
    SSLAnalyserThread*    m_thread;

    enum State {
        Idle,
        Busy,
        Cancelled
    };

    State              m_state;
};

#endif // SSLANALYSER_H

