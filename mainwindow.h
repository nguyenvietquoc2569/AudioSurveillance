#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "engine.h"
#include "levelmeter.h"
#include "settingsdialog.h"
#include "SSLGraph.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    void on_BTSetting_clicked();

    void on_BTStart_clicked();

    void on_BTStop_clicked();

    void bg_collectDataCompleted();
    void bg_gmmInitComplete(Panda_GMM *,qint64 pos,qint64 length);

    void SlotBackgrounddetected();
    void SlotForegrounddetected(qint64,qint64,QString);

    void on_BTSaveBackground_clicked();

private:
    Ui::MainWindow *ui;
    Engine* engine;
    SettingsDialog *m_settingdiag;
    LevelMeter* levelmeter1;
    LevelMeter* levelmeter2;
    LevelMeter* levelmeter3;
    LevelMeter* levelmeter4;
    SSLGraph* sslGraph;

};

#endif // MAINWINDOW_H
