/*  HoughLine
    Copyright (C) 2019 Patrick Molenaar (pr_molenaar@hotmail.com)

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#include "houghline.h"

#include <qmath.h>
#include <QPointF>

/**
 * Initialises the hough line
 */
HoughLine::HoughLine(double theta, double r, int width, int height)
{
    this->theta = theta;
    this->r = r;

    // During processing h_h is doubled so that -ve r values
    int houghHeight = (int) (qSqrt(2) * qMax(height, width)) / 2;

    // Find edge points and vote in array
    float centerX = width / 2;
    float centerY = height / 2;

    // Draw edges in output array
    double sinTheta = qSin(theta);
    double cosTheta = qCos(theta);

    double x1=0.0;
    double y1=0.0;
    double x2=0.0;
    double y2=0.0;
    if (theta < M_PI * 0.25 || theta > M_PI * 0.75) {
        // Store coordinates of vertical-ish lines
        y2 = height - 1.0;
        x1 = ((((r - houghHeight) - ((y1 - centerY) * sinTheta)) / cosTheta) + centerX);
        x2 = ((((r - houghHeight) - ((y2 - centerY) * sinTheta)) / cosTheta) + centerX);
    } else {
        // Store coordinates of horizontal-ish lines
        x2 = width - 1.0;
        y1 = ((((r - houghHeight) - ((x1 - centerX) * cosTheta)) / sinTheta) + centerY);
        y2 = ((((r - houghHeight) - ((x2 - centerX) * cosTheta)) / sinTheta) + centerY);
    }
    beginPoint.setX(x1);
    beginPoint.setY(y1);
    endPoint.setX(x2);
    endPoint.setY(y2);
}

QPointF HoughLine::getBeginPoint() const {
    return beginPoint;
}

QPointF HoughLine::getEndPoint() const {
    return endPoint;
}

/**
 * Sources for intersection and distance calculations came from
 *     http://paulbourke.net/geometry/pointlineplane/
 * Also check https://doc.qt.io/archives/qt-4.8/qlinef.html for more line methods!
 */
HoughLine::IntersectResult HoughLine::Intersect(const HoughLine& other_line, QPointF& intersection)
{
    double denom = ((other_line.getEndPoint().y() - other_line.getBeginPoint().y()) * (endPoint.x() - beginPoint.x())) -
            ((other_line.getEndPoint().x() - other_line.getBeginPoint().x()) * (endPoint.y() - beginPoint.y()));

    double nume_a = ((other_line.getEndPoint().x() - other_line.getBeginPoint().x()) * (beginPoint.y() - other_line.getBeginPoint().y())) -
            ((other_line.getEndPoint().y() - other_line.getBeginPoint().y()) * (beginPoint.x() - other_line.getBeginPoint().x()));

    double nume_b = ((endPoint.x() - beginPoint.x()) * (beginPoint.y() - other_line.getBeginPoint().y())) -
            ((endPoint.y() - beginPoint.y()) * (beginPoint.x() - other_line.getBeginPoint().x()));

    if (denom == 0.0f)
    {
        if(nume_a == 0.0f && nume_b == 0.0f)
        {
            return COINCIDENT;
        }
        return PARALLEL;
    }

    double ua = nume_a / denom;
    double ub = nume_b / denom;

    if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
    {
        // Get the intersection point.
        intersection.setX(qreal(beginPoint.x() + ua * (endPoint.x() - beginPoint.x())));
        intersection.setY(qreal(beginPoint.y() + ua * (endPoint.y() - beginPoint.y())));
        return INTERESECTING;
    }

    return NOT_INTERESECTING;
}

double HoughLine::Magnitude(const QPointF& point1, const QPointF& point2)
{
    QPointF vector = point2 - point1;
//    vector.setX(point2.x - point1.x);
//    vector.setY(point2.y - point1.y);
    return qSqrt(vector.x() * vector.x() + vector.y() * vector.y());
}

bool HoughLine::DistancePointLine(const QPointF& point, QPointF& intersection, float& distance)
{
    double lineMag = Magnitude(endPoint, beginPoint);

    double U = qreal((((point.x() - beginPoint.x()) * (endPoint.x() - beginPoint.x())) +
            ((point.y() - beginPoint.y()) * (endPoint.y() - beginPoint.y()))) /
            (lineMag * lineMag));

    if (U < 0.0 || U > 1.0) {
        return false; // closest point does not fall within the line segment
    }

    intersection.setX(beginPoint.x() + U * (endPoint.x() - beginPoint.x()));
    intersection.setY(beginPoint.y() + U * (endPoint.y() - beginPoint.y()));

    distance = Magnitude(point, intersection);

    return true;
}
