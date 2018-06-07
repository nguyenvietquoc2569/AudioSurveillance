
#ifndef SPECTRUMANALYSER_H
#define SPECTRUMANALYSER_H
#define CHECKED_CONNECT(source, signal, receiver, slot) \
    if (!connect(source, signal, receiver, slot)) \
        qt_assert_x(Q_FUNC_INFO, "CHECKED_CONNECT failed", __FILE__, __LINE__);


#include <QByteArray>
#include <QObject>
#include <QVector>
#include "mfccLib.h"

#ifdef DUMP_SPECTRUMANALYSER
#include <QDir>
#include <QFile>
#include <QTextStream>
#endif


#ifndef DISABLE_FFT
#endif

QT_FORWARD_DECLARE_CLASS(QAudioFormat)
QT_FORWARD_DECLARE_CLASS(QThread)

class FFTRealWrapper;

class SpectrumAnalyserThreadPrivate;

/**
 * Implementation of the spectrum analysis which can be run in a
 * separate thread.
 */
class SpectrumAnalyserThread : public QObject
{
    Q_OBJECT

public:
    SpectrumAnalyserThread(QObject *parent,Panda_mfccParameter *para);
    ~SpectrumAnalyserThread();
    Panda_mfccParameter *mf_format;

public slots:

    void calculateSpectrum(const QByteArray &buffer,
                           qint64 pos,
                           qint64 length);

signals:
    void calculationComplete(double *CC,qint64,qint64);



private:

    QThread*                                    m_thread;

};

/**
 * Class which performs frequency spectrum analysis on a window of
 * audio samples, provided to it by the Engine.
 */
class SpectrumAnalyser : public QObject
{
    Q_OBJECT
public:
    SpectrumAnalyser(QObject *parent = 0,Panda_mfccParameter *para=0);
    ~SpectrumAnalyser();
    Panda_mfccParameter *mf_format;
public:
    void calculate(const QByteArray &buffer, qint64 pos, qint64 lenth);
    bool isReady() const;
    void cancelCalculation();
signals:
    void spectrumChanged(double *cc,qint64,qint64);
private slots:
    void calculationComplete(double * cc,qint64,qint64);
private:
    void calculateWindow();
private:
    SpectrumAnalyserThread*    m_thread;

    enum State {
        Idle,
        Busy,
        Cancelled
    };

    State              m_state;
};

#endif // SPECTRUMANALYSER_H

