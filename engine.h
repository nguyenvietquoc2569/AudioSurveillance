#ifndef ENGINE_H
#define ENGINE_H

#define CHECKED_CONNECT(source, signal, receiver, slot) \
    if (!connect(source, signal, receiver, slot)) \
        qt_assert_x(Q_FUNC_INFO, "CHECKED_CONNECT failed", __FILE__, __LINE__);

#include <QObject>
#include <QAudio>
#include <QAudioInput>
#include <QAudioFormat>
#include "mfccLib.h"
#include "spectrumanalyser.h"
#include "pbufferfeature.h"
#include "GMManalyser.h"
#include "backgrounddetection.h"
#include "SSLanalyser.h"


class Engine : public QObject
{
    Q_OBJECT
public:
    explicit Engine(QObject *parent = 0);
    bool Initialize();
    void startRecording();
    void stopRecording();
    bool setFormat(QAudioFormat);
    bool setAudiodevice(const QAudioDeviceInfo &audiodevice1,const QAudioDeviceInfo &audiodevice2);

    const QList<QAudioDeviceInfo> &availableAudioInputDevices() const
                                    { return m_availableAudioInputDevices; }
    QAudioFormat        m_format;

    const QList<QAudioDeviceInfo> m_availableAudioInputDevices;
    QAudioDeviceInfo    m_audioInputDevice_1;
    QAudioInput*        m_audioInput_1;
    QIODevice*          m_audioInputIODevice_1;
    qint64              m_recordPosition_1;

    int                 m_levelBufferLength_1;
    qreal               m_rmsLevel_1;
    qreal               m_peakLevel_1;

    int                 m_levelBufferLength_2;
    qreal               m_rmsLevel_2;
    qreal               m_peakLevel_2;

    int                 m_levelBufferLength_3;
    qreal               m_rmsLevel_3;
    qreal               m_peakLevel_3;

    int                 m_levelBufferLength_4;
    qreal               m_rmsLevel_4;
    qreal               m_peakLevel_4;


    QByteArray          m_buffer_1;
    qint64              m_bufferPosition_1;
    qint64              m_bufferLength_1;
    qint64              m_dataLength_1;

    QAudioDeviceInfo    m_audioInputDevice_2;
    QAudioInput*        m_audioInput_2;
    QIODevice*          m_audioInputIODevice_2;
    qint64              m_recordPosition_2;

    QByteArray          m_buffer_2;
    qint64              m_bufferPosition_2;
    qint64              m_bufferLength_2;
    qint64              m_dataLength_2;

    QByteArray          c_chanel_1;
    qint64              c_pos_1;
    qint64              c_chanellength1;
    qint64              c_datalength1;
    QByteArray          c_chanel_2;
    qint64              c_pos_2;
    qint64              c_chanellength2;
    qint64              c_datalength2;
    QByteArray          c_chanel_3;
    qint64              c_pos_3;
    qint64              c_chanellength3;
    qint64              c_datalength3;
    QByteArray          c_chanel_4;
    qint64              c_pos_4;
    qint64              c_chanellength4;
    qint64              c_datalength4;

    qint64              mf_length;
    qint64              mf_shiftLength;
    Panda_mfccParameter* mf_format;
    qint64              mf_pos;
    int                 mf_atchannel;
    qreal               mf_max;

    qint64              bf_numberOfFeatureForInitBG;
    qint64              bf_numberOfFeatureForFrame;
    qint64              bf_numberOfFeatureForShiftFrame;

    PBufferFeature       bf_background;
    PBufferFeature       bf_feature;

    BackgroundDetection SS_background;

    int                 m_count;

    bool                Synch_yet;
    bool                Synch_au1ready;
    bool                Synch_au2ready;
    qint64              Synch_num;



signals:
    void recordPositionChanged(qint64 m_recordPosition_1,qint64 m_recordPosition_2);
    void dataLengthChanged(qint64,qint64);
    void BG_SignalCollectdataComplete();
    void BG_GMMSignalInitComplete(Panda_GMM *GMM,qint64 pos,qint64 length);
    void levelChanged_1(qreal rmsLevel, qreal peakLevel, int numSamples);
    void levelChanged_2(qreal rmsLevel, qreal peakLevel, int numSamples);
    void levelChanged_3(qreal rmsLevel, qreal peakLevel, int numSamples);
    void levelChanged_4(qreal rmsLevel, qreal peakLevel, int numSamples);

    void SignalBackgroundDetected();
    void SignalForegroundDetected(qint64 pos,qint64 length,QString classname);
    void SignalSSLChanged(double angel1,double angel2,QString);
    
public slots:
    void audioDataReady_1();
    void audioDataReady_2();
    void audioNotify();
    void MFCCChanged(double *,qint64,qint64);
    void SSLChanged(double angel1,double angel2,qint64,qint64);
    void BG_CollectdataComplete(double** arr,int num,qint64 pos,qint64 length);
    void BG_GMMInitComplete(Panda_GMM *GMM,qint64 pos,qint64 length);

    void SlotBackgroundDetected(qint64 pos,qint64 length);
    void SlotForegroundDetected(qint64 pos,qint64 length,QString classname);


    
private:
    SpectrumAnalyser    *m_spectrumAnalyser;
    GMMAnalyser         *m_GMMAnalyser;
    SSLAnalyser         *m_SSLAnalyser;
    void setRecordPosition(qint64 position_1,qint64 position_2, bool forceEmit);
    qreal calculateLevel_1(qint64,qint64);
    qreal calculateLevel_2(qint64,qint64);
    qreal calculateLevel_3(qint64,qint64);
    qreal calculateLevel_4(qint64,qint64);
    //qreal calculateLevel_1(qint64,qint64);
    double* calculateMFCC1(qint64,qint64);
    double* calculateMFCC2(qint64,qint64);
    double* calculateMFCC3(qint64,qint64);
    double* calculateMFCC4(qint64,qint64);

    void setLevel_1(qreal rmsLevel, qreal peakLevel, int numSamples);
    void setLevel_2(qreal rmsLevel, qreal peakLevel, int numSamples);
    void setLevel_3(qreal rmsLevel, qreal peakLevel, int numSamples);
    void setLevel_4(qreal rmsLevel, qreal peakLevel, int numSamples);
    void setdefaultFormat();

    enum{
        CollectBackground,
        TrainBackGround,
        Working
    } state;

};

#endif // ENGINE_H
