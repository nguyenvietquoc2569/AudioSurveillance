#ifndef BUFFERFEATURE_H
#define BUFFERFEATURE_H

#include <QObject>
#include <QList>
#include <QQueue>

class bufferFeature : public QList<double*>
{
    Q_OBJECT
public:
    explicit bufferFeature(QObject *parent = 0);
    void psetpmaxlength(qint64);
    qint64 pgetpmaxlength();
    void psetpminlength(qint64 length);
    qint64 pgetpminlength();
    void perasep();
    void padd(double*);
    
signals:
    void pfullplength();
    void pminplength();

public slots:

private:
    qint64 pmaxlength;
    qint64 pminlength;
};

#endif // BUFFERFEATURE_H
