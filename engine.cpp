#include "engine.h"
#include "utils.h"
#include "QApplication"
#include "math.h"
#include "mfccLib.h"
#include "spectrumanalyser.h"

//-----------------------------------------------------------------------------
// Constants
//-----------------------------------------------------------------------------
const qreal SSLIntervalMili         = 400;

const qint64 BufferDurationUs       = 100 * 1000000;
const int    NotifyIntervalMs       = 100;
const int    LevelWindowUs          = 0.1 * 1000000;
const qreal  MFCCIntervalMs         = 30;
const qreal  MFCCShiftIntervalMs    = 20;
const long   f_low                  = 100;
const long   f_sampleRate           = 44100;
const long   f_high                 = 20000;
const int    mf_c_M                 = 20;
const int    mf_c_N                 = 22;
const int    mf_c_L                 = 30;
const double mf_c_alpha             = 0.97f;


const int    bf_timeBackgroundInitInMili = 2000;
const int    bf_timeforframe             = 500;
const int    bf_timeforShiftframe        = 400;

const qreal  bg_threshold        = 30;
//------
const double baseline            = 0.13f;
const double soundspeed            = 350;
const double alpha            = 0;
//
const QString FolderClassifier = "D:/ModelForSoundClassification";
const QString FolderBackgroundModel = "D:/Modelforbackground";

Engine::Engine(QObject *parent)
    :    QObject(parent)
    ,   m_audioInputDevice_1(QAudioDeviceInfo::defaultInputDevice())
    ,   m_audioInput_1(0)
    ,   m_audioInputIODevice_1(0)
    ,   m_recordPosition_1(0)
    ,   m_bufferPosition_1(0)
    ,   m_bufferLength_1(0)
    ,   m_dataLength_1(0)
    ,   m_availableAudioInputDevices
            (QAudioDeviceInfo::availableDevices(QAudio::AudioInput))
    ,   m_audioInputDevice_2(QAudioDeviceInfo::defaultInputDevice())
    ,   m_audioInput_2(0)
    ,   m_audioInputIODevice_2(0)
    ,   m_recordPosition_2(0)
    ,   m_bufferPosition_2(0)
    ,   m_bufferLength_2(0)
    ,   m_dataLength_2(0)
    ,   m_count(0)

    ,   m_levelBufferLength_1(0)
    ,   m_rmsLevel_1(0.0)
    ,   m_peakLevel_1(0.0)
    ,   m_levelBufferLength_2(0)    
    ,   m_rmsLevel_2(0.0)    
    ,   m_peakLevel_2(0.0)
    ,   m_levelBufferLength_3(0)
    ,   m_rmsLevel_3(0.0)
    ,   m_peakLevel_3(0.0)
    ,   m_levelBufferLength_4(0)
    ,   m_rmsLevel_4(0.0)    
    ,   m_peakLevel_4(0.0)

    ,   c_pos_1(0)
    ,   c_pos_2(0)
    ,   c_pos_3(0)
    ,   c_pos_4(0)
    ,   c_chanellength1(0)
    ,   c_chanellength2(0)
    ,   c_chanellength3(0)
    ,   c_chanellength4(0)
    ,   c_datalength1(0)
    ,   c_datalength2(0)
    ,   c_datalength3(0)
    ,   c_datalength4(0)
    ,   bf_numberOfFeatureForInitBG(0)
    ,   bf_numberOfFeatureForFrame(0)
    ,   bf_numberOfFeatureForShiftFrame(0)
    ,   bf_background()
    ,   bf_feature()
    ,   state(CollectBackground)
    ,   SS_background()
{

}
bool Engine::Initialize()
{
    state=CollectBackground;
    setdefaultFormat();

    //Configure MFCC
    mf_length   = audioLengthFromTime(m_format,MFCCIntervalMs);
    mf_shiftLength  = audioLengthFromTime(m_format,MFCCShiftIntervalMs);
    mf_format = Panda_Mfcc_CreateParam((long)m_format.sampleRate(),(long)mf_length,(long)mf_shiftLength,"hamming",(long)f_low,(long)f_high,(int)mf_c_M,(int)mf_c_N,(int)mf_c_L,0.97);
    mf_pos = 0;

    m_spectrumAnalyser = new SpectrumAnalyser(this,mf_format);
    m_GMMAnalyser = new GMMAnalyser(this,mf_format,FolderBackgroundModel);
    qDebug() <<"SSL length"<< (1<<FFTSSL::log2(audioLengthFromTime(m_format,MFCCIntervalMs*10)));
    m_SSLAnalyser = new SSLAnalyser(this,baseline,soundspeed,361,alpha,1<<FFTSSL::log2(audioLengthFromTime(m_format,MFCCIntervalMs*10)),audioLengthFromTime(m_format,MFCCShiftIntervalMs*10),m_format.sampleRate());


    SS_background.mf_format=mf_format;
    SS_background.threshold = bg_threshold;
    SS_background.ReadClassifierFromFolder(FolderClassifier);


    CHECKED_CONNECT(m_spectrumAnalyser, SIGNAL(spectrumChanged(double*,qint64,qint64)),
                        this, SLOT(MFCCChanged(double*,qint64,qint64)));
    CHECKED_CONNECT(m_GMMAnalyser, SIGNAL(GMMInitComplete_background(Panda_GMM*,qint64,qint64)),
                    this, SLOT(BG_GMMInitComplete(Panda_GMM *,qint64,qint64)));
    CHECKED_CONNECT(&bf_background, SIGNAL(pfullplength(double**,int,qint64,qint64)),
                    this, SLOT(BG_CollectdataComplete(double**,int,qint64,qint64)));

    CHECKED_CONNECT(this, SIGNAL(BG_GMMSignalInitComplete(Panda_GMM*,qint64,qint64)),
                    &SS_background, SLOT(GMMChange(Panda_GMM*,qint64,qint64)));

    CHECKED_CONNECT(&bf_feature, SIGNAL(pfullplength(double**,int,qint64,qint64)),
                    &SS_background, SLOT(FeaturebufferChange(double**,int,qint64,qint64)));

    CHECKED_CONNECT(&SS_background, SIGNAL(BackgroundDetected(qint64,qint64)),
                    this, SLOT(SlotBackgroundDetected(qint64,qint64)));
    CHECKED_CONNECT(&SS_background, SIGNAL(ForegroundDetected(qint64,qint64,QString)),
                    this, SLOT(SlotForegroundDetected(qint64,qint64,QString)));
    CHECKED_CONNECT(m_SSLAnalyser, SIGNAL(SSLChanged(double,double,qint64,qint64)),
                    this, SLOT(SSLChanged(double,double,qint64,qint64)));


    //Configure buffer for background
    bf_numberOfFeatureForInitBG = NumberOfFeature(m_format,mf_format,bf_timeBackgroundInitInMili);
    bf_numberOfFeatureForFrame = NumberOfFeature(m_format,mf_format,bf_timeforframe);
    bf_numberOfFeatureForShiftFrame = NumberOfFeature(m_format,mf_format,bf_timeforShiftframe);

    bf_background.psetpmaxlength(bf_numberOfFeatureForInitBG);
    bf_background.psetpminlength(bf_numberOfFeatureForInitBG/2);
    bf_feature.psetpmaxlength(bf_numberOfFeatureForFrame);
    bf_feature.psetpminlength(bf_numberOfFeatureForFrame-bf_numberOfFeatureForShiftFrame);
    bf_feature.clear();
    bf_background.clear();


    if (m_audioInputDevice_1.deviceName()==m_audioInputDevice_2.deviceName())
    {
        qDebug()<<"Initialize :: "<<"audio input device have not set";
        return false;
    }
    if (!(m_audioInputDevice_1.isFormatSupported(m_format) & m_audioInputDevice_2.isFormatSupported(m_format)))
    {
        qDebug()<<"Initialize :: "<<"audio input device have not been supported";
        return false;
    }
    delete m_audioInput_1;
    m_audioInput_1 = 0;
    m_audioInputIODevice_1 = 0;

    delete m_audioInput_2;
    m_audioInput_2 = 0;
    m_audioInputIODevice_2 = 0;

    m_recordPosition_1 = 0;
    m_recordPosition_2 = 0;

    m_bufferLength_1 = audioLength(m_format, BufferDurationUs);
    m_buffer_1.resize(m_bufferLength_1);
    m_buffer_1.fill(0);

    m_bufferLength_2 = audioLength(m_format, BufferDurationUs);
    m_buffer_2.resize(m_bufferLength_2);
    m_buffer_2.fill(0);

    c_chanellength1 = audioLength(m_format, BufferDurationUs)/m_format.channelCount();
    c_chanel_1.resize(c_chanellength1);
    c_chanellength2 = audioLength(m_format, BufferDurationUs)/m_format.channelCount();
    c_chanel_2.resize(c_chanellength2);
    c_chanellength3 = audioLength(m_format, BufferDurationUs)/m_format.channelCount();
    c_chanel_3.resize(c_chanellength3);
    c_chanellength4 = audioLength(m_format, BufferDurationUs)/m_format.channelCount();
    c_chanel_4.resize(c_chanellength4);

    c_chanel_1.fill(0);c_datalength1=0;
    c_chanel_2.fill(0);c_datalength2=0;
    c_chanel_3.fill(0);c_datalength3=0;
    c_chanel_4.fill(0);c_datalength4=0;


    m_audioInput_1 = new QAudioInput(m_audioInputDevice_1, m_format, this);
    m_audioInput_1->setNotifyInterval(NotifyIntervalMs);

    m_audioInput_2 = new QAudioInput(m_audioInputDevice_2, m_format, this);
    m_audioInput_2->setNotifyInterval(NotifyIntervalMs);

    qDebug()<<"Initialize :: "<<"succesful"<<c_chanel_2.size()<<c_chanel_1.size();
    return true;
}
bool Engine::setAudiodevice(const QAudioDeviceInfo &audiodevice1,const QAudioDeviceInfo &audiodevice2)
{
    m_audioInputDevice_1 = audiodevice1;
    m_audioInputDevice_2 = audiodevice2;
    qDebug()<<"audioSetdevice :: choose  "<<audiodevice1.deviceName()<<" "<<audiodevice2.deviceName();
    return true;
}
void Engine::startRecording()
{
    if (!m_audioInput_1 || !m_audioInput_1)
    {
        qDebug()<<"StartRecording :: error "<<m_audioInput_1<<" "<<m_audioInput_2;
    }
    m_buffer_1.fill(0);
    m_buffer_2.fill(0);
    setRecordPosition(0,0, true);
//    CHECKED_CONNECT(m_audioInput_1, SIGNAL(stateChanged(QAudio::State)),
//                    this, SLOT(audioStateChanged(QAudio::State)));
//    CHECKED_CONNECT(m_audioInput_2, SIGNAL(stateChanged(QAudio::State)),
//                    this, SLOT(audioStateChanged(QAudio::State)));
    CHECKED_CONNECT(m_audioInput_1, SIGNAL(notify()),
                    this, SLOT(audioNotify()));

    //CHECKED_CONNECT(m_audioInput_2, SIGNAL(notify()),this, SLOT(audioNotify()));

    m_count = 0;
    m_dataLength_1 = 0;
    m_dataLength_2 = 0;
    emit dataLengthChanged(0,0);
    m_audioInputIODevice_1 = m_audioInput_1->start();
    m_audioInputIODevice_2 = m_audioInput_2->start();
    CHECKED_CONNECT(m_audioInputIODevice_1, SIGNAL(readyRead()),
                    this, SLOT(audioDataReady_1()));
    CHECKED_CONNECT(m_audioInputIODevice_2, SIGNAL(readyRead()),
                    this, SLOT(audioDataReady_2()));
    Synch_yet = false;
    Synch_au1ready =false;
    Synch_au2ready =false;
    Synch_num =0 ;
    qDebug()<<"StartRecording :: successful ";
}

void Engine::stopRecording()
{
    if (m_audioInput_1 || m_audioInput_2) {
        m_audioInput_1->stop();
        m_audioInput_2->stop();
        QCoreApplication::instance()->processEvents();
        m_audioInput_1->disconnect();
        m_audioInput_2->disconnect();
    }
    m_audioInputIODevice_1 = 0;
    m_audioInputIODevice_2 = 0;
    Synch_au1ready =false;
    Synch_au2ready = false;
    Synch_num =0;
    Synch_yet =false;
    qDebug()<<"stopped";
}

bool Engine::setFormat(QAudioFormat format)
{

    if (true == m_audioInputDevice_1.isFormatSupported(format))
    {
        m_format = format;
        qDebug()<<"SetFormat :: "<<"successful";
        m_levelBufferLength_1 = audioLength(m_format,LevelWindowUs);
        m_levelBufferLength_2 = audioLength(m_format,LevelWindowUs);
        m_levelBufferLength_3 = audioLength(m_format,LevelWindowUs);
        m_levelBufferLength_4 = audioLength(m_format,LevelWindowUs);
        qDebug()<<"SetFormat :: "<<"set Level Window" << m_levelBufferLength_1<<m_levelBufferLength_2;
        return true;
    }
    qDebug()<<"SetFormat :: "<<"fail";
    return false;
}
//Private
void Engine::setRecordPosition(qint64 position1,qint64 position2, bool forceEmit)
{
    const bool changed = (m_recordPosition_1 != position1 || m_recordPosition_2 != position2);
    m_recordPosition_1 = position1;
    m_recordPosition_2 = position2;
    if (changed || forceEmit)
        emit recordPositionChanged(m_recordPosition_1,m_recordPosition_2);
}
//Slot

void Engine::audioDataReady_1()
{
    Q_ASSERT(0 == m_bufferPosition_1);
    Synch_au1ready = true;

    const qint64 bytesReady = m_audioInput_1->bytesReady();

    if (!Synch_au2ready)
    {
        m_audioInputIODevice_1->read(bytesReady);
        //qDebug()<<"Waiting device 2 , reject " <<bytesReady;
        return;
    }


    const qint64 bytesSpace = m_buffer_1.size() - m_dataLength_1;
    const qint64 bytesToRead = qMin(bytesReady, bytesSpace);

    const qint64 bytesRead = m_audioInputIODevice_1->read(
                                       m_buffer_1.data() + m_dataLength_1,
                                       bytesToRead);

    qint64 step = m_format.channelCount()*m_format.sampleSize()/8;
    for(qint64 j=0;j<bytesRead;j=j+step)
    {
        memcpy(c_chanel_1.data()+c_datalength1,
                 m_buffer_1.data()+m_dataLength_1+j,
               m_format.sampleSize()/8
                    );
        memcpy(c_chanel_2.data()+c_datalength2,
               m_buffer_1.data()+m_dataLength_1+j+m_format.sampleSize()/8,
               m_format.sampleSize()/8
               );
        c_datalength1=c_datalength1+m_format.sampleSize()/8;
        c_datalength2=c_datalength2+m_format.sampleSize()/8;
    }




    if (bytesRead) {
        m_dataLength_1 += bytesRead;
        emit dataLengthChanged(m_dataLength_1,m_dataLength_2);
    }
    if (m_buffer_1.size() == m_dataLength_1 || c_chanel_1.size()==c_datalength1 || c_chanel_2.size()==c_datalength2)
        stopRecording();


}

void Engine::audioDataReady_2()
{
    Q_ASSERT(0 == m_bufferPosition_2);
    Synch_au2ready=true;

    //qDebug()<<"byte doc rang tu 1" <<m_audioInput_1->bytesReady();

    const qint64 bytesReady = m_audioInput_2->bytesReady();
    const qint64 bytesSpace = m_buffer_2.size() - m_dataLength_2;
    const qint64 bytesToRead = qMin(bytesReady, bytesSpace);

    const qint64 bytesRead = m_audioInputIODevice_2->read(
                                       m_buffer_2.data() + m_dataLength_2,
                                       bytesToRead);


    qint64 step = m_format.channelCount()*m_format.sampleSize()/8;
    for(qint64 j=0;j<bytesRead;j=j+step)
    {
        memcpy(c_chanel_3.data()+c_datalength3,
                 m_buffer_2.data()+m_dataLength_2+j,
               m_format.sampleSize()/8
                    );
        memcpy(c_chanel_4.data()+c_datalength4,
               m_buffer_2.data()+m_dataLength_2+j+m_format.sampleSize()/8,
               m_format.sampleSize()/8
               );
        c_datalength3=c_datalength3+m_format.sampleSize()/8;
        c_datalength4=c_datalength4+m_format.sampleSize()/8;
    }

    if (bytesRead) {
        m_dataLength_2 += bytesRead;
        emit dataLengthChanged(m_dataLength_1,m_dataLength_2);
    }

    if (m_buffer_2.size() == m_dataLength_2 || c_chanel_3.size()==c_datalength3 || c_chanel_4.size()==c_datalength4)
        stopRecording();
}

void Engine::audioNotify()
{
    qreal temp;
    const qint64 recordPosition_1 = qMin(m_bufferLength_1, audioLength(m_format, m_audioInput_1->processedUSecs()));
    const qint64 recordPosition_2 = qMin(m_bufferLength_2, audioLength(m_format, m_audioInput_2->processedUSecs()));
    setRecordPosition(recordPosition_1,recordPosition_2,false);

    //const qint64 levelPosition_1 = m_dataLength_1 - m_levelBufferLength_1;
    const qint64 levelPosition_1 = c_datalength1 - m_levelBufferLength_1;
    const qint64 levelPosition_2 = c_datalength2 - m_levelBufferLength_2;
    const qint64 levelPosition_3 = c_datalength3 - m_levelBufferLength_3;
    const qint64 levelPosition_4 = c_datalength4 - m_levelBufferLength_4;

    mf_max=0;
    mf_atchannel = 1;

    if (levelPosition_1 >= mf_max)
    {
        temp = calculateLevel_1(levelPosition_1, m_levelBufferLength_1);
        if (temp > mf_max) { mf_max = temp; mf_atchannel = 1;}
    }
    if (levelPosition_2 >= mf_max)
    {
        temp = calculateLevel_2(levelPosition_2, m_levelBufferLength_2);
        if (temp > mf_max) { mf_max = temp; mf_atchannel = 2;}
    }
    if (levelPosition_3 >= mf_max)
    {
        temp = calculateLevel_3(levelPosition_3, m_levelBufferLength_3);
        if (temp > mf_max) { mf_max = temp; mf_atchannel = 3;}
    }
    if (levelPosition_4 >= mf_max)
    {
        temp = calculateLevel_4(levelPosition_4, m_levelBufferLength_4);
        if (temp > mf_max) { mf_max = temp; mf_atchannel = 4;}
    }
    while (mf_pos + mf_length*m_format.sampleSize()/8 <= c_pos_1 +c_datalength1 &&
           mf_pos + mf_length*m_format.sampleSize()/8 <= c_pos_2 +c_datalength2 &&
           mf_pos + mf_length*m_format.sampleSize()/8 <= c_pos_3 +c_datalength3 &&
           mf_pos + mf_length*m_format.sampleSize()/8 <= c_pos_4 +c_datalength4 )
    {
        if (mf_atchannel==1)
        {
//            double *mfcc = calculateMFCC1(mf_pos,mf_length);
//            qDebug()<< mf_pos << mfcc[0]<< mfcc[1]<< mfcc[2];
            m_spectrumAnalyser->calculate(c_chanel_1,mf_pos,mf_length);
            mf_pos = mf_pos+mf_shiftLength*m_format.sampleSize()/8;

        }
        else if (mf_atchannel==2)
        {
//            double *mfcc = calculateMFCC2(mf_pos,mf_length);
//            qDebug()<< mf_pos << mfcc[0]<< mfcc[1]<< mfcc[2];
            m_spectrumAnalyser->calculate(c_chanel_2,mf_pos,mf_length);
            mf_pos = mf_pos+mf_shiftLength*m_format.sampleSize()/8;
        }
        else if (mf_atchannel==3)
        {
//            double *mfcc = calculateMFCC3(mf_pos,mf_length);
//            qDebug()<< mf_pos << mfcc[0]<< mfcc[1]<< mfcc[2];
            m_spectrumAnalyser->calculate(c_chanel_3,mf_pos,mf_length);
            mf_pos = mf_pos+mf_shiftLength*m_format.sampleSize()/8;
        }
        else if (mf_atchannel==4)
        {
//            double *mfcc = calculateMFCC4(mf_pos,mf_length);
//            qDebug()<< mf_pos << mfcc[0]<< mfcc[1]<< mfcc[2];
            m_spectrumAnalyser->calculate(c_chanel_4,mf_pos,mf_length);
            mf_pos = mf_pos+mf_shiftLength*m_format.sampleSize()/8;
        }
    }

}

void Engine::MFCCChanged(double * mfcc,qint64 pos,qint64 length)
{
    //qDebug() <<pos<<length<< mfcc[0]<< mfcc[1]<< mfcc[2];

    if (state==CollectBackground)
        bf_background.padd(mfcc,pos,length);

    if (state==Working)
        bf_feature.padd(mfcc,pos,length);
}

void Engine::SSLChanged(double angel1, double angel2, qint64 pos, qint64 length)
{
    SignalSSLChanged(angel1,angel2,"unknow");
    qDebug()<<"("<<angel1<<","<<angel2 <<") at "<<pos;
}

void Engine::BG_CollectdataComplete(double** arr,int num,qint64 pos,qint64 length)
{
    qDebug()<<"Engine:: Slot BG_CollectdataComplete : ";
    qDebug()<<"         State change to Training Background ";
    state=TrainBackGround;
    m_GMMAnalyser->calculate_backgournd(bf_background,(bf_background.listpos.first()),(bf_background.listpos.last()+bf_background.listlen.last()) );
    emit BG_SignalCollectdataComplete();

}

void Engine::BG_GMMInitComplete(Panda_GMM *GMM,qint64 pos,qint64 length)
{
    qDebug()<<"Engine:: Slot BG_GMMInitComplete : ";
    qDebug()<<"         State change to Working ";
    emit BG_GMMSignalInitComplete(GMM,pos,length);
    state = Working;
}

void Engine::SlotBackgroundDetected(qint64 pos,qint64 length)
{
    emit SignalBackgroundDetected();
}

void Engine::SlotForegroundDetected(qint64 pos,qint64 length,QString classname)
{
    m_SSLAnalyser->calculate(c_chanel_3,c_chanel_4,c_chanel_1,c_chanel_2,pos,length);
    emit SignalForegroundDetected(pos,length,classname);
}


///////////////////Private
qreal Engine::calculateLevel_1(qint64 position, qint64 length)
{
#ifdef DISABLE_LEVEL
    Q_UNUSED(position)
    Q_UNUSED(length)
#else
    Q_ASSERT(position + length <= c_pos_1 + c_chanellength1);

    qreal peakLevel = 0.0;

    qreal sum = 0.0;
    const char *ptr = c_chanel_1.constData() + position - c_pos_1;
    const char *const end = ptr + length;
    while (ptr < end) {
        const qint16 value = *reinterpret_cast<const qint16*>(ptr);
        const qreal fracValue = pcmToReal(value);
        peakLevel = qMax(peakLevel, fracValue);
        sum += fracValue * fracValue;
        ptr += 2;
    }
    const int numSamples = length / 2;
    qreal rmsLevel = sqrt(sum / numSamples);

    rmsLevel = qMax(qreal(0.0), rmsLevel);
    rmsLevel = qMin(qreal(1.0), rmsLevel);
    setLevel_1(rmsLevel, peakLevel, numSamples);

//    qDebug() << "Engine::calculateLevel 1" << "pos" << position << "len" << length
//                 << "rms" << rmsLevel << "peak" << peakLevel;
    return rmsLevel;
#endif
}
qreal Engine::calculateLevel_2(qint64 position, qint64 length)
{
#ifdef DISABLE_LEVEL
    Q_UNUSED(position)
    Q_UNUSED(length)
#else
    Q_ASSERT(position + length <= c_pos_2 + c_chanellength2);

    qreal peakLevel = 0.0;

    qreal sum = 0.0;
    const char *ptr = c_chanel_2.constData() + position - c_pos_2;
    const char *const end = ptr + length;
    while (ptr < end) {
        const qint16 value = *reinterpret_cast<const qint16*>(ptr);
        const qreal fracValue = pcmToReal(value);
        peakLevel = qMax(peakLevel, fracValue);
        sum += fracValue * fracValue;
        ptr += 2;
    }
    const int numSamples = length / 2;
    qreal rmsLevel = sqrt(sum / numSamples);

    rmsLevel = qMax(qreal(0.0), rmsLevel);
    rmsLevel = qMin(qreal(1.0), rmsLevel);
    setLevel_2(rmsLevel, peakLevel, numSamples);

//    qDebug() << "Engine::calculateLevel 2" << "pos" << position << "len" << length
//                 << "rms" << rmsLevel << "peak" << peakLevel;
    return rmsLevel;
#endif
}
qreal Engine::calculateLevel_3(qint64 position, qint64 length)
{
#ifdef DISABLE_LEVEL
    Q_UNUSED(position)
    Q_UNUSED(length)
#else
    Q_ASSERT(position + length <= c_pos_3 + c_chanellength3);

    qreal peakLevel = 0.0;

    qreal sum = 0.0;
    const char *ptr = c_chanel_3.constData() + position - c_pos_3;
    const char *const end = ptr + length;
    while (ptr < end) {
        const qint16 value = *reinterpret_cast<const qint16*>(ptr);
        const qreal fracValue = pcmToReal(value);
        peakLevel = qMax(peakLevel, fracValue);
        sum += fracValue * fracValue;
        ptr += 2;
    }
    const int numSamples = length / 2;
    qreal rmsLevel = sqrt(sum / numSamples);

    rmsLevel = qMax(qreal(0.0), rmsLevel);
    rmsLevel = qMin(qreal(1.0), rmsLevel);
    setLevel_3(rmsLevel, peakLevel, numSamples);

//    qDebug() << "Engine::calculateLevel 3" << "pos" << position << "len" << length
//                 << "rms" << rmsLevel << "peak" << peakLevel;
    return rmsLevel;
#endif
}
qreal Engine::calculateLevel_4(qint64 position, qint64 length)
{
#ifdef DISABLE_LEVEL
    Q_UNUSED(position)
    Q_UNUSED(length)
#else
    Q_ASSERT(position + length <= c_pos_4 + c_chanellength4);

    qreal peakLevel = 0.0;

    qreal sum = 0.0;
    const char *ptr = c_chanel_4.constData() + position - c_pos_4;
    const char *const end = ptr + length;
    while (ptr < end) {
        const qint16 value = *reinterpret_cast<const qint16*>(ptr);
        const qreal fracValue = pcmToReal(value);
        peakLevel = qMax(peakLevel, fracValue);
        sum += fracValue * fracValue;
        ptr += 2;
    }
    const int numSamples = length / 2;
    qreal rmsLevel = sqrt(sum / numSamples);

    rmsLevel = qMax(qreal(0.0), rmsLevel);
    rmsLevel = qMin(qreal(1.0), rmsLevel);
    setLevel_4(rmsLevel, peakLevel, numSamples);

//    qDebug() << "Engine::calculateLevel 4" << "pos" << position << "len" << length
//                 << "rms" << rmsLevel << "peak" << peakLevel;
    return rmsLevel;
#endif
}

double* Engine::calculateMFCC1(qint64 position ,qint64 length)
{
    return PandaSampleToMfcc(mf_format,(short*)(c_chanel_1.data()+position),length);
}
double* Engine::calculateMFCC2(qint64 position ,qint64 length)
{
    return PandaSampleToMfcc(mf_format,(short*)(c_chanel_2.data()+position),length);
}
double* Engine::calculateMFCC3(qint64 position ,qint64 length)
{
    return PandaSampleToMfcc(mf_format,(short*)(c_chanel_3.data()+position),length);
}
double* Engine::calculateMFCC4(qint64 position ,qint64 length)
{
    return PandaSampleToMfcc(mf_format,(short*)(c_chanel_4.data()+position),length);
}


void Engine::setLevel_1(qreal rmsLevel, qreal peakLevel, int numSamples)
{
    m_rmsLevel_1 = rmsLevel;
    m_peakLevel_1 = peakLevel;
    emit levelChanged_1(m_rmsLevel_1, m_peakLevel_1, numSamples);
}
void Engine::setLevel_2(qreal rmsLevel, qreal peakLevel, int numSamples)
{
    m_rmsLevel_2 = rmsLevel;
    m_peakLevel_2 = peakLevel;
    emit levelChanged_2(m_rmsLevel_2, m_peakLevel_2, numSamples);
}

void Engine::setLevel_3(qreal rmsLevel, qreal peakLevel, int numSamples)
{
    m_rmsLevel_3 = rmsLevel;
    m_peakLevel_3 = peakLevel;
    emit levelChanged_3(m_rmsLevel_3, m_peakLevel_3, numSamples);
}
void Engine::setLevel_4(qreal rmsLevel, qreal peakLevel, int numSamples)
{
    m_rmsLevel_4 = rmsLevel;
    m_peakLevel_4 = peakLevel;
    emit levelChanged_4(m_rmsLevel_4, m_peakLevel_4, numSamples);
}
void Engine::setdefaultFormat()
{
    QAudioFormat format;
    format.setByteOrder(QAudioFormat::LittleEndian);
    format.setCodec("audio/pcm");
    format.setSampleSize(16);
    format.setSampleType(QAudioFormat::SignedInt);
    format.setSampleRate(f_sampleRate );
    format.setChannelCount(2);
    setFormat(format);
}
