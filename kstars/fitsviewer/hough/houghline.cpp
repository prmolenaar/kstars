/*  HoughLine
    Copyright (C) 2019 Patrick Molenaar (pr_molenaar@hotmail.com)

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#include "houghline.h"

#include <qmath.h>

/**
 * Initialises the hough line
 */
HoughLine::HoughLine(double theta, double r, int width, int height, int score)
    : QLineF()
{
    this->theta = theta;
    this->r = r;
    this->score = score;

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
    setP1(QPointF(x1, y1));
    setP2(QPointF(x2, y2));
}

int HoughLine::getScore() const
{
    return score;
}

double HoughLine::getR() const
{
    return r;
}

double HoughLine::getTheta() const
{
    return theta;
}

void HoughLine::setTheta(const double theta)
{
    this->theta = theta;
}

bool HoughLine::compareByScore(const HoughLine *line1,const HoughLine *line2)
{
    return (line1->getScore() < line2->getScore());
}

bool HoughLine::compareByTheta(const HoughLine *line1,const HoughLine *line2)
{
    return (line1->getTheta() < line2->getTheta());
}

/**
 * Sources for intersection and distance calculations came from
 *     http://paulbourke.net/geometry/pointlineplane/
 * Also check https://doc.qt.io/archives/qt-4.8/qlinef.html for more line methods
 */
HoughLine::IntersectResult HoughLine::Intersect(const HoughLine& other_line, QPointF& intersection)
{
    double denom = ((other_line.p2().y() - other_line.p1().y()) * (p2().x() - p1().x())) -
            ((other_line.p2().x() - other_line.p1().x()) * (p2().y() - p1().y()));

    double nume_a = ((other_line.p2().x() - other_line.p1().x()) * (p1().y() - other_line.p1().y())) -
            ((other_line.p2().y() - other_line.p1().y()) * (p1().x() - other_line.p1().x()));

    double nume_b = ((p2().x() - p1().x()) * (p1().y() - other_line.p1().y())) -
            ((p2().y() - p1().y()) * (p1().x() - other_line.p1().x()));

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
        intersection.setX(qreal(p1().x() + ua * (p2().x() - p1().x())));
        intersection.setY(qreal(p1().y() + ua * (p2().y() - p1().y())));
        return INTERESECTING;
    }

    return NOT_INTERESECTING;
}

double HoughLine::Magnitude(const QPointF& point1, const QPointF& point2)
{
    QPointF vector = point2 - point1;
    return qSqrt(vector.x() * vector.x() + vector.y() * vector.y());
}

bool HoughLine::DistancePointLine(const QPointF& point, QPointF& intersection, double& distance)
{
    double lineMag = length();

    double U = qreal((((point.x() - p1().x()) * (p2().x() - p1().x())) +
            ((point.y() - p1().y()) * (p2().y() - p1().y()))) /
            (lineMag * lineMag));

    if (U < 0.0 || U > 1.0) {
        return false; // closest point does not fall within the line segment
    }

    intersection.setX(p1().x() + U * (p2().x() - p1().x()));
    intersection.setY(p1().y() + U * (p2().y() - p1().y()));

    distance = Magnitude(point, intersection);

    return true;
}
