/*  FITS Hough Transform
    Copyright (C) 2019 Patrick Molenaar (pr_molenaar@hotmail.com)

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#ifndef HOUGHTRANSFORM_H_
#define HOUGHTRANSFORM_H_

#include "houghline.h"

#include <math.h>
#include <QVector>
#include <QRgb>

/**
 * @class HoughTransform
 * Hough transform implementation
 * Based on the java implementation found on http://vase.essex.ac.uk/software/HoughTransform
 * @author Patrick Molenaar
 * @version 1.0
 */
class HoughTransform
{
  public:
    template <typename T>
    HoughTransform(int width, int height, QVector<T> &houghArray);

    template <typename T>
    void initialise(QVector<T> &houghArray);

    template <typename T>
    void addPoints(QVector<T> &image, QVector<T> &houghArray);

    template <typename T>
    void addPoint(int x, int y, QVector<T> &houghArray);

    template <typename T>
    void getLines(int threshold, QVector<T> &houghArray, QVector<HoughLine*> &lines);

    template <typename T>
    int getHighestValue(QVector<T> &houghArray);

    template <typename T>
    QVector<QRgb> getHoughArrayImage(QVector<T> &houghArray);

  private:
    // The size of the neighbourhood in which to search for other local maxima
    int neighbourhoodSize = 4;

    // How many discrete values of theta shall we check?
    int maxTheta = 180;

    // Using maxTheta, work out the step
    double thetaStep = M_PI / maxTheta;

    // the width and height of the image
    int width;
    int height;

    // the hough array (size: width x height)
//    QVector<T> &houghArray;

    // the coordinates of the centre of the image
    float centerX, centerY;

    // the height of the hough array
    int houghHeight;

    // double the hough height (allows for negative numbers)
    int doubleHeight;

    // the number of points that have been added
    int numPoints;

    // cache of values of sin and cos for different theta values. Has a significant performance improvement.
    QVector<double> sinCache;
    QVector<double> cosCache;
};

#endif
