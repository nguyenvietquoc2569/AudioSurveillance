#include "bufferfeature.h"
#include <QObject>
#include <QList>
#include <QQueue>
bufferFeature::bufferFeature(QObject *parent) :
    QList<double*>()
{
    pmaxlength=1;
    pminlength=0;
}

void bufferFeature::psetpmaxlength(qint64 length)
{
    pmaxlength = length;
}

qint64 bufferFeature::pgetpmaxlength()
{
    return pmaxlength;
}

void bufferFeature::psetpminlength(qint64 length)
{
    pminlength = length;
}

qint64 bufferFeature::pgetpminlength()
{
    return pminlength;
}

void bufferFeature::perasep()
{
    while (count()>pminlength)
    {
        double *d = first();
        free(d);
        removeFirst();
    }
    emit pminplength();
}

void bufferFeature::padd(double *feature)
{
    push_back(feature);
    if (count()>= pmaxlength)
        emit pfullplength();
}
