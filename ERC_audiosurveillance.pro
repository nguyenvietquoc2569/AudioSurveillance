#-------------------------------------------------
#
# Project created by QtCreator 2013-04-27T16:22:03
#
#-------------------------------------------------

QT       += core gui multimedia

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ERC_audiosurveillance
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    utils.cpp \
    settingsdialog.cpp \
    levelmeter.cpp \
    engine.cpp \
    complex_simple.cpp \
    dct.cpp \
    Panda_GMM_Util.cpp \
    Wav_tool.cpp \
    clust_invert.cpp \
    invert.cpp \
    alloc_util.cpp \
    mfccLib.cpp \
    Panda_GMM.cpp \
    spectrumanalyser.cpp \
    pbufferfeature.cpp \
    GMManalyser.cpp \
    backgrounddetection.cpp \
    fft.cpp \
    fftssl.cpp \
    SSLanalyser.cpp \
    SSLGraph.cpp \
    wavfile.cpp

HEADERS  += mainwindow.h \
    utils.h \
    settingsdialog.h \
    levelmeter.h \
    engine.h \
    dct.h \
    alloc_util.h \
    invert.h \
    WavTool.h \
    mfccLib.h \
    Panda_GMM_Util.h \
    Panda_GMM.h \
    complex_simple.h \
    spectrumanalyser.h \
    pbufferfeature.h \
    GMManalyser.h \
    backgrounddetection.h \
    fft.h \
    fftssl.h \
    SSLanalyser.h \
    SSLGraph.h \
    wavfile.h

FORMS    += mainwindow.ui
