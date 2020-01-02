/*  HoughTransform
    Copyright (C) 2019 Patrick Molenaar (pr_molenaar@hotmail.com)

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#include "houghtransform.h"

#include <qmath.h>
#include <fits_debug.h>

/**
 * Initialises the hough transform. The dimensions of the input image are needed
 * in order to initialise the hough array.
 *
 * @param width  The width of the input image
 * @param height The height of the input image
 */
template <typename T>
HoughTransform::HoughTransform(int width, int height, QVector<T> &houghArray)
{
    this->width = width;
    this->height = height;
    initialise(houghArray);
}

/**
 * Initialises the hough array. Called by the constructor so you don't need to call it
 * yourself, however you can use it to reset the transform if you want to plug in another
 * image (although that image must have the same width and height)
 */
template <typename T>
void HoughTransform::initialise(QVector<T> &houghArray)
{
    // Calculate the maximum height the hough array needs to have
    houghHeight = (int) (qSqrt(2) * qMax(height, width)) / 2;

    // Double the height of the hough array to cope with negative r values
    doubleHeight = 2 * houghHeight;

    // Create the hough array (resize to MAX_THETA * doubleHeight and fill with 0)
    houghArray.fill((T)0, MAX_THETA * doubleHeight);

    // Find edge points and vote in array
    centerX = width / 2;
    centerY = height / 2;

    // Count how many points there are
    numPoints = 0;

    // cache the values of sin and cos for faster processing
    sinCache.resize(MAX_THETA);
    cosCache.resize(MAX_THETA);
    for (int t = 0; t < MAX_THETA; t++)
    {
        double realTheta = t * thetaStep;
        sinCache[t] = qSin(realTheta);
        cosCache[t] = qCos(realTheta);
    }
}

/**
 * Adds points from an image. The image is assumed to be greyscale black and white, so all pixels that are
 * not black are counted as edges.
 */
template <typename T>
void HoughTransform::addPoints(const QVector<T> &image, QVector<T> &houghArray) {

    // Now find edge points and update the hough array
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            // Find non-black pixels
            int index = y * width + x;
            if (image[index] != 0)
            {
                addPoint(x, y, houghArray);
            }
        }
    }
}

/**
 * Adds a single point to the hough transform. You can use this method directly
 * if your data isn't represented as a buffered image.
 */
template <typename T>
void HoughTransform::addPoint(int x, int y, QVector<T> &houghArray) {

    // Go through each value of theta
    for (int t = 0; t < MAX_THETA; t++)
    {
        //Work out the r values for each theta step
        int r = (int) (((x - centerX) * cosCache[t]) + ((y - centerY) * sinCache[t]));

        // this copes with negative values of r
        r += houghHeight;

        if (r < 0 || r >= doubleHeight)
            continue;

        // Increment the hough array
        int index = r * MAX_THETA + t; // houghArray[t][r]
        if (index >= 0 && index < houghArray.size())
        {
            houghArray[index]++;
        }
        else
        {
            qCDebug(KSTARS_FITS) << "Error: Index out of range " << index << " should be smaller then " << houghArray.size();
        }
    }
    numPoints++;
}

/**
 * Once points have been added in some way this method extracts the lines and returns them as a Vector
 * of HoughLine objects, which can be used to draw on the
 *
 * @param percentageThreshold The percentage threshold above which lines are determined from the hough array
 */
template <typename T>
void HoughTransform::getLines(int threshold, QVector<T> &houghArray, QVector<HoughLine*> &lines)
{
    // Only proceed if the hough array is not empty
    if (numPoints == 0)
    {
        return;
    }

    // Search for local peaks above threshold to draw
    for (int t = 0; t < MAX_THETA; t++)
    {
        bool loop = false;
        for (int r = NEIGHBOURHOOD_SIZE; r < doubleHeight - NEIGHBOURHOOD_SIZE && !loop; r++)
        {
            // Only consider points above threshold
            int index = r * MAX_THETA + t; // houghArray[t][r]
            if (houghArray[index] > threshold)
            {
                int peak = houghArray[index];

                // Check that this peak is indeed the local maxima
                for (int dx = -NEIGHBOURHOOD_SIZE; dx <= NEIGHBOURHOOD_SIZE && !loop; dx++)
                {
                    for (int dy = -NEIGHBOURHOOD_SIZE; dy <= NEIGHBOURHOOD_SIZE && !loop; dy++)
                    {
                        int dt = t + dx;
                        int dr = r + dy;
                        if (dt < 0) dt = dt + MAX_THETA;
                        else if (dt >= MAX_THETA) dt = dt - MAX_THETA;
                        int dIndex = dr * MAX_THETA + dt; // houghArray[t][r]
                        if (houghArray[dIndex] > peak) {
                            // found a bigger point nearby, skip
                            loop = true;
                        }
                    }
                }

                if (!loop)
                {
                    // calculate the true value of theta
                    double theta = t * thetaStep;

                    // add the line to the vector
                    HoughLine *pHoughLine = new HoughLine(theta, r, width, height, peak);
                    lines.append(pHoughLine);
                }
            }
        }
    }
}

void HoughTransform::getSortedTopThreeLines(QVector<HoughLine*> &houghLines, QVector<HoughLine*> &top3Lines)
{
    // Sort houghLines by score (highest scores are clearest lines)
    // TODO PRM: Use of sort compare methods see: https://www.off-soft.net/en/develop/qt/qtb1.html
    qSort(houghLines.begin(),houghLines.end(), HoughLine::compareByScore);

    // Get top three lines (these should represent the three lines matching the bahtinov mask lines
    top3Lines = houghLines.mid(0, 3);
    // Verify the angle of these lines with regard to the bahtinov mask angle, correct the angle if necessary
    HoughLine* lineR = top3Lines[0];
    HoughLine* lineG = top3Lines[1];
    HoughLine* lineB = top3Lines[2];
    double thetaR = lineR->getTheta();
    double thetaG = lineG->getTheta();
    double thetaB = lineB->getTheta();
    // Calculate angle between each line
    double bahtinovMaskAngle = qDegreesToRadians(20.0 + 5.0); // bahtinov mask angle plus 5 degree margin
    double dGR = thetaR - thetaG;
    double dBG = thetaB - thetaG;
    double dBR = thetaB - thetaR;
    if (dGR > bahtinovMaskAngle && dBR > bahtinovMaskAngle) {
        // lineR has theta that is 180 degrees rotated
        thetaR -= M_PI;
        // update theta
        lineR->setTheta(thetaR);
    }
    if (dBR > bahtinovMaskAngle && dBG > bahtinovMaskAngle) {
        // lineB has theta that is 180 degrees rotated
        thetaB -= M_PI;
        // update theta
        lineB->setTheta(thetaB);
    }
    if (dGR > bahtinovMaskAngle && dBG > bahtinovMaskAngle) {
        // lineG has theta that is 180 degrees rotated
        thetaG -= M_PI;
        // update theta
        lineG->setTheta(thetaG);
    }
    // Now sort top3lines array according to calculated new angles
    qSort(top3Lines.begin(),top3Lines.end(), HoughLine::compareByTheta);
}

/**
 * Gets the highest value in the hough array
 */
template <typename T>
int HoughTransform::getHighestValue(QVector<T> &houghArray)
{
    int max = 0;
    for (int t = 0; t < MAX_THETA; t++)
    {
        for (int r = 0; r < doubleHeight; r++)
        {
            int index = r * MAX_THETA + t; // houghArray[t][r]
            if (houghArray[index] > max)
            {
                max = houghArray[index];
            }
        }
    }
    return max;
}

/**
 * Gets the hough array as an image, in case you want to have a look at it.
 */
template <typename T>
QVector<QRgb> HoughTransform::getHoughArrayImage(QVector<T> &houghArray)
{
    int max = getHighestValue(houghArray);
    QVector<QRgb> image = new QVector<QRgb>(MAX_THETA * doubleHeight);
    for (int t = 0; t < MAX_THETA; t++)
    {
        for (int r = 0; r < doubleHeight; r++)
        {
            int index = r * MAX_THETA + t; // houghArray[t][r]
            double value = 255 * ((double) houghArray[index]) / max;
            int v = 255 - (int) value;
            QRgb c = qRgb(v, v, v);
            image[index] = c;
        }
    }
    return image;
}
