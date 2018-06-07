
#include "SSLanalyser.h"
#include "utils.h"
#include <qmath.h>
#include <qmetatype.h>
#include <QAudioFormat>
#include <QThread>
#include "mfccLib.h"
#include "fftSSL.h"
const qint64 framesize = 2734;
const double PI  =3.141592653589793238462;
SSLAnalyserThread::SSLAnalyserThread(QObject *parent, double baseline,double soundspeed, int nDelays, double beta, qint64 FFTSize,qint64 shiftFrame,qint64 SampleRate)
    :   QObject(parent)
    ,   m_thread(new QThread(this))

{

    // moveToThread() cannot be called on a QObject with a parent
    setParent(0);
    double mMaxITDsecs = (double)baseline / (double)soundspeed;
    mFFTSize = FFTSize;
    mFFTSize_2 = FFTSize>>1;
    mShiftFrame = shiftFrame;
    mBeta = beta;
    mNDelays = nDelays;
    msoundspeed = soundspeed;
    mSamplerate = SampleRate;
    mFFt = new FFTSSL(mFFTSize);
    // Create window mask
    mWindow = new double[mFFTSize];
    double factor = 2.0 * PI / (mFFTSize - 1);
    double f = 0.0;
    for (int i = 0; i < mFFTSize; ++i)
    {
           mWindow[i] = 0.54f - 0.46f * (double) qCos(f);
                f += factor;
    }
    //create complex phase delay
    // init complex phase delays
    mDelayFactors = new double[2 * (mNDelays * mFFTSize_2)];
    for (int i = 0; i < mNDelays; ++i) {
        // delay value in seconds
        double delay = mMaxITDsecs / 2.0 * qSin(i / (mNDelays - 1.0) * PI - PI / 2.0);
        for (int j = 0; j < mFFTSize_2; ++j) {
            // frequency (in Hz) corresponding to FFT index j
            double frequency = j * mSamplerate/ ((double) mFFTSize);
            // phase angle corresponding to delay value
            double phaseAngle = -2.0 * PI * frequency * delay;
            // store real part
            mDelayFactors[(i * mFFTSize_2 + j) << 1] = (double) qCos(phaseAngle);
            // store imaginary part
            mDelayFactors[((i * mFFTSize_2 + j) << 1) + 1] = (double) qSin(phaseAngle);
        }
    }
    mCoinMap = new double[(mNDelays * mFFTSize_2)];
    //initialize AZithm
    initAz();
    moveToThread(m_thread);
    m_thread->start();
    m_thread->setPriority(QThread::HighPriority);

}

SSLAnalyserThread::~SSLAnalyserThread()
{
}

double SSLAnalyserThread::CalculateForeachPart(const QByteArray &buffer1,const QByteArray &buffer2,qint64 pos,qint64 length)
{
    double ketqua;
    length = length / 2;

    memset(mCoinMap,0,sizeof(double)*(mNDelays * mFFTSize_2));

    short* left = (short*) (buffer1.data()+pos);
    short* right = (short*) (buffer2.data()+pos);
    int i = 0;

    double *dleft = new double[mFFTSize];
    double *dright = new double[mFFTSize];

    //qDebug()<<"begin to FFT";
    qreal MaxOfEnergy = -1.0f;
    while (i+mFFTSize-1<=length)
    {
        qreal temp = CalcEnergy(left,i,mFFTSize);
        if (temp>MaxOfEnergy)
        {
            MaxOfEnergy = temp;
            memset(dleft,0,sizeof(double)*mFFTSize);
            memset(dright,0,sizeof(double)*mFFTSize);

            int z=0;
            for(int j=i;j<i+mFFTSize;j++)
            {
                dleft[z]=pcmToReal((qint16)left[j])*mWindow[z];
                dright[z]=pcmToReal((qint16)right[j])*mWindow[z];
                z++;
            }
            //qDebug()<<dleft[0]<<dright[0];
            double *fftleft = mFFt->compute(dleft,mFFTSize);
            double *fftright = mFFt->compute(dright,mFFTSize);


            calculatemCoinc(fftleft,fftright);
            free(fftleft);
            free(fftright);
        }
        i=i+mShiftFrame;
    }

    double* result = new double[mNDelays];

    double max=-1.79769313486231e308;
    double min= 1.79769313486232e308;

    for (i = 0; i < mNDelays; ++i)
    {
        for (int j = 0; j < mFFTSize_2; ++j) {
            result[mNDelays - 1 - i] += mCoinMap[i * mFFTSize_2 + j];
        }
        if (result[mNDelays - 1 - i] >= max) {
            max = result[mNDelays - 1 - i];
        } else if (result[mNDelays - 1 - i] < min) {
            min = result[mNDelays - 1 - i];
        }
    }

    double normfac = 1.0f / (max - min);
    for (i = 0; i < mNDelays; ++i)
        result[i] = (max - result[i]) * normfac;

    int i1 = mNDelays - 1, i2 = 1;
    for (i = 0; i < mNDelays; ++i)
    {
        if ((result[i] > 0.99f) && (result[i] > result[i1])&& (result[i] > result[i2]))
                    ketqua=mAzVals[i];
        i1 = (i1 + 1) % mNDelays;
        i2 = (i2 + 1) % mNDelays;
    }

    free(dleft);
    free(dright);
    free(result);
    return ketqua;
}

void SSLAnalyserThread::calculatemCoinc(double *leftFFT, double *rightFFT)
{
    for (int i = 0; i < mNDelays; ++i)
    {
        for (int j = 100; j < mFFTSize_2/4; ++j)
        {
            double ReL = leftFFT[j << 1], ImL = leftFFT[(j << 1) + 1];
            double ReR = rightFFT[j << 1], ImR = rightFFT[(j << 1) + 1];
            double ReD = mDelayFactors[(i * mFFTSize_2 + j) << 1];
            double ImD = mDelayFactors[((i * mFFTSize_2 + j) << 1) + 1];
            double delRe = (ReL - ReR) * ReD - (ImL + ImR) * ImD;
            double delIm = (ImL - ImR) * ReD + (ReL + ReR) * ImD;
            mCoinMap[i * mFFTSize_2 + j] = mCoinMap[i * mFFTSize_2 + j]
                    * mBeta + (double) qSqrt(delRe * delRe + delIm * delIm);
        }
    }
}

void SSLAnalyserThread::initAz()
{
    mAzVals = new double[mNDelays];
    //mAzVals =(double*)malloc(sizeof(double)*mNDelays);
    for (int i = 0; i < mNDelays; ++i)
    {
        mAzVals[mNDelays - 1 - i] = 90.0f - i / (mNDelays - 1.0f) * 180.0f;
    }
}

void SSLAnalyserThread::calculateSSL(const QByteArray &buffer,
                                     const QByteArray &buffer1,
                                     const QByteArray &buffer2,
                                     const QByteArray &buffer3,
                                  qint64 pos,
                                  qint64 length)
{
    double re1=CalculateForeachPart(buffer,buffer1,pos,length);
    double re2=CalculateForeachPart(buffer2,buffer3,pos,length);
    emit calculationComplete(re1,re2,pos,length);
}


//=============================================================================
// SSLAnalyser
//=============================================================================

SSLAnalyser::SSLAnalyser(QObject *parent, double baseline,double soundspeed, int nDelays, double beta, qint64 FFTSize,qint64 shiftFrame,qint64 SampleRate)
    :   QObject(parent)    
    ,   m_state(Idle)
#ifdef DUMP_SSLANALYSER
    ,   m_count(0)

#endif
{
    m_thread= new SSLAnalyserThread(this,baseline,soundspeed,nDelays,beta,FFTSize,shiftFrame,SampleRate);
    CHECKED_CONNECT(m_thread, SIGNAL(calculationComplete(double,double,qint64,qint64)),
                    this, SLOT(calculationComplete(double,double,qint64,qint64)));
}

SSLAnalyser::~SSLAnalyser()
{ 
}

//-----------------------------------------------------------------------------
// Public functions
//-----------------------------------------------------------------------------


void SSLAnalyser::calculate(const QByteArray &buffer,
                            const QByteArray &buffer1,
                            const QByteArray &buffer2,
                            const QByteArray &buffer3,
                         qint64 pos,
                         qint64 length)
{
    // QThread::currentThread is marked 'for internal use only', but
    // we're only using it for debug output here, so it's probably OK :)
//    SSLANALYSER_DEBUG << "SSLAnalyser::calculate"
//                           << QThread::currentThread()
//                           << "state" << m_state;

    //if (isReady()) {
        m_state = Busy;

        // Invoke SSLAnalyserThread::calculateSSL using QMetaObject.  If
        // m_thread is in a different thread from the current thread, the
        // calculation will be done in the child thread.
        // Once the calculation is finished, a calculationChanged signal will be
        // emitted by m_thread.
        const bool b = QMetaObject::invokeMethod(m_thread, "calculateSSL",
                                                 Qt::AutoConnection,
                                                 Q_ARG(QByteArray, buffer),
                                                 Q_ARG(QByteArray, buffer1),
                                                 Q_ARG(QByteArray, buffer2),
                                                 Q_ARG(QByteArray, buffer3),
                                                 Q_ARG(qint64, pos),
                                                 Q_ARG(qint64, length));
        Q_ASSERT(b);
        Q_UNUSED(b) // suppress warnings in release builds
    //}
}

bool SSLAnalyser::isReady() const
{
    return (Idle == m_state);
}

void SSLAnalyser::cancelCalculation()
{
    if (Busy == m_state)
        m_state = Cancelled;
}


//-----------------------------------------------------------------------------
// Private slots
//-----------------------------------------------------------------------------

void SSLAnalyser::calculationComplete(double angel1,double angel2,qint64 pos,qint64 length)
{
    //Q_ASSERT(Idle != m_state);
    //if (Busy == m_state)
        emit SSLChanged(angel1,angel2,pos,length);
    //m_state = Idle;
}
