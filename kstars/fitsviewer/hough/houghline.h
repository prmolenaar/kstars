/*  HoughLine
    Copyright (C) 2019 Patrick Molenaar (pr_molenaar@hotmail.com)

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#ifndef HOUGHLINE_H_
#define HOUGHLINE_H_

#include <QtCore/QObject>
#include <QPointF>

/**
 * @class HoughLine
 * Line representation for HoughTransform
 * Based on the java implementation found on http://vase.essex.ac.uk/software/HoughTransform
 * @author Patrick Molenaar
 * @version 1.0
 */
class HoughLine : public QObject
{
    Q_OBJECT

  public:

    enum IntersectResult {
        PARALLEL,
        COINCIDENT,
        NOT_INTERESECTING,
        INTERESECTING
    };

    HoughLine(double theta, double r, int width, int height);
    IntersectResult Intersect(const HoughLine& other_line, QPointF& intersection);
    bool DistancePointLine(const QPointF& point, QPointF& intersection, float& distance);
    QPointF getBeginPoint() const;
    QPointF getEndPoint() const;

  private:
    double Magnitude(const QPointF& point1, const QPointF& point2);

    double theta;
    double r;
    QPointF beginPoint;
    QPointF endPoint;
};

#endif
