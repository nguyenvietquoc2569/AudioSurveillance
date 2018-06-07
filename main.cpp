#include "mainwindow.h"
#include <QApplication>
#include <QDebug>
#include "pbufferfeature.h"
#include <QDir>
#include "fftSSL.h"
int main(int argc, char *argv[])
{
    qRegisterMetaType<PBufferFeature>("PBufferFeature");
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
