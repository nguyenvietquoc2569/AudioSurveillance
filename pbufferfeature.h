#ifndef PBUFFERFEATURE_H
#define PBUFFERFEATURE_H

#include <QObject>

class PBufferFeature : public QObject
{
    Q_OBJECT
public:
    explicit PBufferFeature(QObject *parent = 0);
    explicit PBufferFeature(const PBufferFeature &buffer);

    void psetpmaxlength(qint64);
    qint64 pgetpmaxlength();
    void psetpminlength(qint64 length);
    qint64 pgetpminlength();
    void clear();
    qint64 count();

    void perasep();
    void padd(double*,qint64 pos,qint64 length);

    double** ToArray();
    QList<double*> list;
    QList<qint64> listpos;
    QList<qint64> listlen;

signals:
    void pfullplength(double **arr,int num,qint64 pos,qint64 length);
    void pminplength();
    
public slots:
    
private:
    qint64 pmaxlength;
    qint64 pminlength;

};
//Q_DECLARE_METATYPE(PBufferFeature)

#endif // PBUFFERFEATURE_H

