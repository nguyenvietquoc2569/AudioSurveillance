/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "SSLGraph.h"

#include <math.h>

#include <QPainter>
#include <QTimer>
#include <QDebug>


// Constants
const int RedrawInterval = 100; // ms
const qreal PeakDecayRate = 0.001;
const int PeakHoldLevelDuration = 2000; // ms


SSLGraph::SSLGraph(QWidget *parent)
    :   QWidget(parent),
      m_redrawTimer(new QTimer(this))

{
    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
    setMinimumWidth(30);

    connect(m_redrawTimer, SIGNAL(timeout()), this, SLOT(redrawTimerExpired()));
    m_redrawTimer->start(RedrawInterval);
    Classname="";
    reset();

}

SSLGraph::~SSLGraph()
{
}

void SSLGraph::reset()
{
    Azimuth = 0;
    Atitude = 0;
    update();
}

void SSLGraph::levelChanged(qreal Azimuth, qreal Atitude, QString Label)
{
    this->Azimuth = Azimuth;
    this->Atitude = Atitude;
    update();
}

void SSLGraph::redrawTimerExpired()
{

    update();
}

void SSLGraph::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event)

    QPainter painter(this);
    painter.fillRect(rect(), Qt::black);

    QPen pen(Qt::white);
    pen.setWidth(1);
    painter.setPen(pen);

    //qDebug()<<rect().center().x()<<0<<rect().center().x()<<rect().bottom();
    painter.drawLine(rect().center().x(),0,rect().center().x(),rect().bottom());
    painter.drawLine(0,rect().center().y(),rect().right(),rect().center().y());
    painter.drawArc(rect().center().x()-10,rect().center().y()-10,20,20,0*16,360*16);

    int delta = 40;
    QString str;

    for (int i = 10;i<=90;i=i+10)
    {
        //int y = rect().center().y;
        //int x = rect().center().x;

        //y - i/10*delta;
        painter.drawLine(rect().center().x()-2,
                         rect().center().y() - i/10*delta,
                         rect().center().x()+2,
                         rect().center().y() - i/10*delta);
        painter.drawLine(rect().center().x()-2,
                         rect().center().y() + i/10*delta,
                         rect().center().x()+2,
                         rect().center().y() + i/10*delta);
        painter.drawLine(rect().center().x()+ i/10*delta,
                         rect().center().y() - 2 ,
                         rect().center().x()+ i/10*delta,
                         rect().center().y() + 2 );
        painter.drawLine(rect().center().x()- i/10*delta,
                         rect().center().y() - 2 ,
                         rect().center().x()- i/10*delta,
                         rect().center().y() + 2 );
        painter.drawText(rect().center().x()- i/10*delta -6,
                         rect().center().y() + 12 ,
                         QString::number(-i));
        painter.drawText(rect().center().x() + i/10*delta -6,
                         rect().center().y() + 12 ,
                         QString::number(i));

        painter.drawText(rect().center().x()-15 ,
                         rect().center().y() - i/10*delta+6 ,
                         QString::number(i));
        painter.drawText(rect().center().x()-18 ,
                         rect().center().y() + i/10*delta+5 ,
                         QString::number(-i));
    }

    pen.setStyle(Qt::DashLine);
    painter.setPen(pen);

    painter.drawLine((int)(rect().center().x()+Atitude/10.0f*delta),
                     rect().center().y(),
                     (int)(rect().center().x()+Atitude/10.0f*delta),
                     (int)(rect().center().y()+Azimuth/10.0f*delta));
    painter.drawLine(rect().center().x(),
                     (int)(rect().center().y()+Azimuth/10.0f*delta),
                     (int)(rect().center().x()+Atitude/10.0f*delta),
                     (int)(rect().center().y()+Azimuth/10.0f*delta));





    QPen slpen(Qt::green);
    slpen.setWidth(3);
    painter.setPen(slpen);

    painter.drawLine(rect().center().x(),
                     rect().center().y(),
                     (int)(rect().center().x()+Atitude/10.0f*delta),
                     (int)(rect().center().y()+Azimuth/10.0f*delta));

    slpen.setColor(Qt::red);
    painter.setPen(slpen);
    painter.drawArc((int)(rect().center().x()+Atitude/10.0f*delta)-3,
                    (int)(rect().center().y()+Azimuth/10.0f*delta)-3,
                          6,
                          6,
                          0*16,
                          360*16);

    QFont font=painter.font() ;
    font.setPointSize ( 20 );
    painter.setFont(font);
    painter.drawText(QPoint((int)(rect().center().x()+Atitude/10.0f*delta),(int)(rect().center().y()+Azimuth/10.0f*delta)),Classname);
    //qDebug()<<"Ve : "<<Atitude<<Azimuth;
}
