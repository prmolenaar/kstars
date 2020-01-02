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

    HoughLine(double theta, double r, int width, int height, int score, QObject *parent = 0);
    IntersectResult Intersect(const HoughLine& other_line, QPointF& intersection);
    bool DistancePointLine(const QPointF& point, QPointF& intersection, double& distance);
    QPointF getBeginPoint() const;
    QPointF getEndPoint() const;
    int getScore() const;
    double getR() const;
    double getTheta() const;
    void setTheta(const double theta);

	// TODO PRM: Check if these methods are necessary!
    HoughLine(const HoughLine &other)
        : QObject(other.parent())
    {
        *this = other;
    }
    HoughLine& operator=(const HoughLine &other)
    {
        theta = other.getTheta();
        r = other.getR();
        score = other.getScore();
        beginPoint = other.getBeginPoint();
        endPoint = other.getEndPoint();
        return *this;
    }
    bool operator<(const HoughLine &other) const
    {
        return (score < other.getScore());
    }

    static bool compareByScore(const HoughLine *line1,const HoughLine *line2);
    static bool compareByTheta(const HoughLine *line1,const HoughLine *line2);

    void printHoughLine()
    {
        printf("Houghline: [s:%d, r:%.2f, theta:%.2f, p1:%.2f,%.2f, p1:%.2f,%.2f]\r\n",
               score, r, theta, beginPoint.x(), beginPoint.y(), endPoint.x(), endPoint.y());
    }

  private:
    double Magnitude(const QPointF& point1, const QPointF& point2);

    int score;
    double theta;
    double r;
    QPointF beginPoint;
    QPointF endPoint;
};

#endif
