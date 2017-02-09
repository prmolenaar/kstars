/*  Ekos Alignment View
 *  Child of AlignView with few additions necessary for Alignment functions

    Copyright (C) 2017 Jasem Mutlaq <mutlaqja@ikarustech.com>

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#include <QtConcurrent>

#include "alignview.h"

#define ZOOM_DEFAULT	100.0
#define ZOOM_MIN        10
#define ZOOM_MAX        400

AlignView::AlignView(QWidget *parent, FITSMode mode, FITSScale filter) : FITSView(parent, mode, filter)
{

}

AlignView::~AlignView()
{

}

void AlignView::drawOverlay(QPainter *painter)
{
    painter->setOpacity(0.4);
    FITSView::drawOverlay(painter);
    painter->setOpacity(1);

    if (RACircle.isNull() == false)
        drawCircle(painter);

    if (correctionLine.isNull() == false)
        drawLine(painter);
}

bool AlignView::createWCSFile(const QString & newWCSFile, double orientation, double ra, double dec, double pixscale)
{
    return imageData->createWCSFile(newWCSFile, orientation, ra, dec, pixscale);
}

void AlignView::setCorrectionParams(QLineF line)
{
    correctionLine   = line;
    correctionCenter = line.p2();

    markerCrosshair.setX(correctionCenter.x());
    markerCrosshair.setY(correctionCenter.y());

    updateFrame();
}

void AlignView::setCorrectionOffset(QPointF newOffset)
{
    if (newOffset.isNull() == false)
    {
        double offsetX = newOffset.x() - imageData->getWidth()/2;
        double offsetY = newOffset.y() - imageData->getHeight()/2;

        correctionOffset.setX(offsetX);
        correctionOffset.setY(offsetY);

        markerCrosshair.setX(correctionCenter.x() + correctionOffset.x());
        markerCrosshair.setY(correctionCenter.y() + correctionOffset.y());
    }
    // Clear points
    else
    {
        correctionOffset = newOffset;
        markerCrosshair  = newOffset;
    }

    updateFrame();
}

void AlignView::drawLine(QPainter *painter)
{
    painter->setPen(QPen( Qt::magenta , 2));
    painter->setBrush( Qt::NoBrush );
    double zoomFactor = (currentZoom / ZOOM_DEFAULT);

    int offsetX = 0, offsetY=0;

    if (correctionOffset.isNull() == false)
    {
        offsetX = correctionOffset.x();
        offsetY = correctionOffset.y();
    }

    double x1 = (correctionLine.p1().x() + offsetX) * zoomFactor;
    double y1 = (correctionLine.p1().y() + offsetY) * zoomFactor;

    double x2 = (correctionLine.p2().x() + offsetX) * zoomFactor;
    double y2 = (correctionLine.p2().y() + offsetY) * zoomFactor;

    QLineF zoomedLine(x1, y1, x2, y2);

    painter->drawLine(zoomedLine);
}

void AlignView::drawCircle(QPainter *painter)
{
    QPen pen(Qt::green);
    pen.setWidth(2);
    pen.setStyle(Qt::DashLine);
    painter->setPen(pen);
    painter->setBrush( Qt::NoBrush);
    double zoomFactor = (currentZoom / ZOOM_DEFAULT);

    QPointF center(RACircle.x()*zoomFactor, RACircle.y() * zoomFactor);
    double r = RACircle.z() * zoomFactor;

    painter->drawEllipse(center, r/25.0, r/25.0);
    painter->drawText((center.x()+5)*zoomFactor,(center.y()+5)*zoomFactor, i18n("RA Axis"));
    painter->drawEllipse(center, r, r);

}


void AlignView::setRACircle(const QVector3D &value)
{
    RACircle = value;
    updateFrame();
}

