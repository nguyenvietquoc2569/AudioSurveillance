#include "pbufferfeature.h"

PBufferFeature::PBufferFeature(QObject *parent) :
    QObject(parent)
  , list()
  , listpos()
  , listlen()
{
    pmaxlength = 1;
    pminlength = 0;
}

PBufferFeature::PBufferFeature(const PBufferFeature &buffer)
{
    this->list=buffer.list;
    this->listpos=buffer.listpos;
    this->listlen=buffer.listlen;
    this->pmaxlength= buffer.pmaxlength;
    this->pminlength= buffer.pminlength;
}

void PBufferFeature::psetpmaxlength(qint64 length)
{
    pmaxlength = length;
}

qint64 PBufferFeature::pgetpmaxlength()
{
    return pmaxlength;
}

void PBufferFeature::psetpminlength(qint64 length)
{
    pminlength = length;
}

qint64 PBufferFeature::pgetpminlength()
{
    return pminlength;
}

void PBufferFeature::clear()
{
    list.clear();
    listlen.clear();
    listpos.clear();
}

qint64 PBufferFeature::count()
{
    return list.count();
}

void PBufferFeature::perasep()
{
    while (list.count()>pminlength)
    {
        //double **d = list.first();
        //free(*d);
        list.removeFirst();
        listlen.removeFirst();
        listpos.removeFirst();
    }
    emit pminplength();
}

void PBufferFeature::padd(double *feature,qint64 pos,qint64 length)
{
    list.push_back(feature);
    listpos.push_back(pos);
    listlen.push_back(length);

    if (list.count()>= pmaxlength)
    {
        QList<double*>::iterator i;
        int dem =0;
        double **arr = (double**)malloc(count()*sizeof(double*));
        for(i=list.begin();i!=list.end();++i)
        {
            arr[dem]=*i;
            dem++;
        }

        emit pfullplength(arr,dem,(listpos[0]),pos + length -(listpos[0]));
        this->perasep();
    }

}

double **PBufferFeature::ToArray()
{
    QList<double*>::iterator i;
    double** re = (double**)malloc(list.count()*sizeof(double*));
    int dem=0;
    for(i=list.begin();i!=list.end();++i)
    {
        re[dem]=*i;
        dem++;
    }
    return re;
}
