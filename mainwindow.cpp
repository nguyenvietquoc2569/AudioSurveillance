#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QStyle>
#include <QMenu>
#include <QFileDialog>
#include <QTimerEvent>
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    levelmeter1(new LevelMeter(this)),
    levelmeter2(new LevelMeter(this)),
    levelmeter3(new LevelMeter(this)),
    levelmeter4(new LevelMeter(this)),
    sslGraph(new SSLGraph(this)),
    engine(new Engine(this)),
    m_settingdiag(new SettingsDialog(
                      engine->availableAudioInputDevices(),
                      engine->availableAudioInputDevices(),
                      this)),
    ui(new Ui::MainWindow)
{

    ui->setupUi(this);

    ui->BTSaveBackground->setVisible(false);
    ui->LBNotification->setVisible(false);
    ui->TBThreashold->setVisible(false);

    levelmeter1->setGeometry(50,50,15,100);
    levelmeter1->reset();
    levelmeter2->setGeometry(75,50,15,100);
    levelmeter2->reset();
    levelmeter3->setGeometry(100,50,15,100);
    levelmeter3->reset();
    levelmeter4->setGeometry(125,50,15,100);
    levelmeter4->reset();
    sslGraph->setGeometry(175,50,800,800);
    sslGraph->reset();

    CHECKED_CONNECT(engine, SIGNAL(levelChanged_1(qreal,qreal,int)),
            levelmeter1, SLOT(levelChanged(qreal, qreal, int)));
    CHECKED_CONNECT(engine, SIGNAL(levelChanged_2(qreal,qreal,int)),
            levelmeter2, SLOT(levelChanged(qreal, qreal, int)));
    CHECKED_CONNECT(engine, SIGNAL(levelChanged_3(qreal,qreal,int)),
            levelmeter3, SLOT(levelChanged(qreal, qreal, int)));
    CHECKED_CONNECT(engine, SIGNAL(levelChanged_4(qreal,qreal,int)),
            levelmeter4, SLOT(levelChanged(qreal, qreal, int)));
    CHECKED_CONNECT(engine, SIGNAL(BG_SignalCollectdataComplete()),
            this, SLOT(bg_collectDataCompleted()));
    CHECKED_CONNECT(engine, SIGNAL(BG_GMMSignalInitComplete(Panda_GMM*,qint64,qint64)),
                    this, SLOT(bg_gmmInitComplete(Panda_GMM*,qint64,qint64)));
    CHECKED_CONNECT(engine, SIGNAL(SignalBackgroundDetected()),
                    this, SLOT(SlotBackgrounddetected()));
    CHECKED_CONNECT(engine, SIGNAL(SignalForegroundDetected(qint64,qint64,QString)),
                    this, SLOT(SlotForegrounddetected(qint64,qint64,QString)));
    CHECKED_CONNECT(engine, SIGNAL(SignalSSLChanged(double,double,QString)),
                    sslGraph, SLOT(levelChanged(qreal,qreal,QString)));



    QVBoxLayout *windowLayout = new QVBoxLayout(this);
    QScopedPointer<QHBoxLayout> analysisLayout(new QHBoxLayout);

    analysisLayout->addWidget(levelmeter1);
    analysisLayout->addWidget(levelmeter2);
    windowLayout->addLayout(analysisLayout.data());
    analysisLayout.take();
    setLayout(windowLayout);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_BTSetting_clicked()
{
    m_settingdiag->exec();
    if (m_settingdiag->result()==1)
    {
        qDebug()<< "on_BtOpenSetting_clicked :: Clicked Ok";
        engine->setAudiodevice(m_settingdiag->input1Device(),m_settingdiag->input2Device());
        return;
    }
    qDebug()<< "on_BtOpenSetting_clicked :: Clicked cancel";
}

void MainWindow::on_BTStart_clicked()
{
    if (engine->Initialize()) engine->startRecording();
}

void MainWindow::on_BTStop_clicked()
{
    engine->stopRecording();
}

void MainWindow::bg_collectDataCompleted()
{
    ui->LBNotification->setText("Ok finish collect");
}

void MainWindow::bg_gmmInitComplete(Panda_GMM *GMM, qint64 pos, qint64 length)
{
    ui->LBNotification->setText("Ok finish training");
}

void MainWindow::SlotBackgrounddetected()
{
    ui->LBBackground->setStyleSheet("QLabel { color : white; }");
    ui->LBBackground->setText("Background");
}
void MainWindow::SlotForegrounddetected(qint64 pos,qint64 length,QString classname)
{
    ui->LBBackground->setStyleSheet("QLabel { color : red; }");
    ui->LBBackground->setText("Foreground");
    sslGraph->Classname=classname;
}

void MainWindow::on_BTSaveBackground_clicked()
{
    engine->SS_background.threshold = ui->TBThreashold->text().toDouble();
}
