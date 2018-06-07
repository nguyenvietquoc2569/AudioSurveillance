
#include "spectrumanalyser.h"
#include "utils.h"
#include <qmath.h>
#include <qmetatype.h>
#include <QAudioFormat>
#include <QThread>
#include "mfccLib.h"
SpectrumAnalyserThread::SpectrumAnalyserThread(QObject *parent,Panda_mfccParameter *para)
    :   QObject(parent)
    ,   m_thread(new QThread(this))

{
    // moveToThread() cannot be called on a QObject with a parent
    mf_format = para;
    setParent(0);
    m_thread->setPriority(QThread::HighPriority);
    moveToThread(m_thread);
    m_thread->start();
}

SpectrumAnalyserThread::~SpectrumAnalyserThread()
{
}

void SpectrumAnalyserThread::calculateSpectrum(const QByteArray &buffer,
                                                qint64 pos,qint64 length)
{
//    qDebug() << "SpectrumAnalyserThread::calculate"
//                           << QThread::currentThread();

    double * CC;
    CC = PandaSampleToMfcc(mf_format,(short*)(buffer.data()+pos),length);
    emit calculationComplete(CC,pos,length);
}


//=============================================================================
// SpectrumAnalyser
//=============================================================================

SpectrumAnalyser::SpectrumAnalyser(QObject *parent,Panda_mfccParameter *para)
    :   QObject(parent)
    ,   m_thread(new SpectrumAnalyserThread(this,para))
    ,   m_state(Idle)
#ifdef DUMP_SPECTRUMANALYSER
    ,   m_count(0)
#endif
{
    CHECKED_CONNECT(m_thread, SIGNAL(calculationComplete(double *,qint64,qint64)),
                    this, SLOT(calculationComplete(double *,qint64,qint64)));
}

SpectrumAnalyser::~SpectrumAnalyser()
{

}

//-----------------------------------------------------------------------------
// Public functions
//-----------------------------------------------------------------------------


void SpectrumAnalyser::calculate(const QByteArray &buffer,
                         qint64 pos,
                         qint64 length)
{
    // QThread::currentThread is marked 'for internal use only', but
    // we're only using it for debug output here, so it's probably OK :)
//    SPECTRUMANALYSER_DEBUG << "SpectrumAnalyser::calculate"
//                           << QThread::currentThread()
//                           << "state" << m_state;

    //if (isReady()) {
        m_state = Busy;

        // Invoke SpectrumAnalyserThread::calculateSpectrum using QMetaObject.  If
        // m_thread is in a different thread from the current thread, the
        // calculation will be done in the child thread.
        // Once the calculation is finished, a calculationChanged signal will be
        // emitted by m_thread.
        const bool b = QMetaObject::invokeMethod(m_thread, "calculateSpectrum",
                                  Qt::AutoConnection,
                                  Q_ARG(QByteArray, buffer),
                                  Q_ARG(qint64, pos),
                                  Q_ARG(qint64, length));
        Q_ASSERT(b);
        Q_UNUSED(b) // suppress warnings in release builds
    //}
}

bool SpectrumAnalyser::isReady() const
{
    return (Idle == m_state);
}

void SpectrumAnalyser::cancelCalculation()
{
    if (Busy == m_state)
        m_state = Cancelled;
}


//-----------------------------------------------------------------------------
// Private slots
//-----------------------------------------------------------------------------

void SpectrumAnalyser::calculationComplete(double *cc,qint64 pos,qint64 length)
{
    //Q_ASSERT(Idle != m_state);
    //if (Busy == m_state)
        emit spectrumChanged(cc,pos,length);
    //m_state = Idle;
}
