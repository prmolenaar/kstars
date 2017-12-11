/*  Ekos PHD2 Handler
    Copyright (C) 2016 Jasem Mutlaq <mutlaqja@ikarustech.com>

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
*/

#include "phd2.h"

#include "Options.h"
#include "kspaths.h"
#include "kstars.h"
#include "fitsio.h"
#include "ekos/ekosmanager.h"

#include <KMessageBox>
#include <QImage>

#include <QJsonDocument>
#include <QJsonObject>
#include <QtNetwork/QNetworkReply>

#include <ekos_guide_debug.h>

#define MAX_SET_CONNECTED_RETRIES   3

namespace Ekos
{
PHD2::PHD2()
{
    tcpSocket = new QTcpSocket(this);

    connect(tcpSocket, SIGNAL(readyRead()), this, SLOT(readPHD2()));
    connect(tcpSocket, SIGNAL(error(QAbstractSocket::SocketError)), this,
            SLOT(displayError(QAbstractSocket::SocketError)));

    events["Version"]                 = Version;
    events["LockPositionSet"]         = LockPositionSet;
    events["CalibrationComplete"]     = CalibrationComplete;
    events["StarSelected"]            = StarSelected;
    events["StartGuiding"]            = StartGuiding;
    events["Paused"]                  = Paused;
    events["StartCalibration"]        = StartCalibration;
    events["AppState"]                = AppState;
    events["CalibrationFailed"]       = CalibrationFailed;
    events["CalibrationDataFlipped"]  = CalibrationDataFlipped;
    events["LoopingExposures"]        = LoopingExposures;
    events["LoopingExposuresStopped"] = LoopingExposuresStopped;
    events["SettleBegin"]             = SettleBegin;
    events["Settling"]                = Settling;
    events["SettleDone"]              = SettleDone;
    events["StarLost"]                = StarLost;
    events["GuidingStopped"]          = GuidingStopped;
    events["Resumed"]                 = Resumed;
    events["GuideStep"]               = GuideStep;
    events["GuidingDithered"]         = GuidingDithered;
    events["LockPositionLost"]        = LockPositionLost;
    events["Alert"]                   = Alert;
}

PHD2::~PHD2()
{
}

bool PHD2::Connect()
{
    if (connection == DISCONNECTED)
    {
        connection = CONNECTING;
        tcpSocket->connectToHost(Options::pHD2Host(), Options::pHD2Port());
    }
    // Already connected, let's connect equipment
    else
        setEquipmentConnected(true);

    return true;
}

bool PHD2::Disconnect()
{
    if (connection == EQUIPMENT_CONNECTED)
        setEquipmentConnected(false);

    connection = DISCONNECTED;
    tcpSocket->disconnectFromHost();

    emit newStatus(GUIDE_DISCONNECTED);

    return true;
}

void PHD2::displayError(QAbstractSocket::SocketError socketError)
{
    switch (socketError)
    {
        case QAbstractSocket::RemoteHostClosedError:
            break;
        case QAbstractSocket::HostNotFoundError:
            emit newLog(i18n("The host was not found. Please check the host name and port settings in Guide options."));
            emit newStatus(GUIDE_DISCONNECTED);
            break;
        case QAbstractSocket::ConnectionRefusedError:
            emit newLog(i18n("The connection was refused by the peer. Make sure the PHD2 is running, and check that "
                             "the host name and port settings are correct."));
            emit newStatus(GUIDE_DISCONNECTED);
            break;
        default:
            emit newLog(i18n("The following error occurred: %1.", tcpSocket->errorString()));
    }

    connection = DISCONNECTED;
}

void PHD2::readPHD2()
{
    QTextStream stream(tcpSocket);

    QJsonParseError qjsonError;

    while (stream.atEnd() == false)
    {
        QString rawString = stream.readLine();

        if (rawString.isEmpty())
            continue;

        QJsonDocument jdoc = QJsonDocument::fromJson(rawString.toLatin1(), &qjsonError);

        if (qjsonError.error != QJsonParseError::NoError)
        {
            emit newLog(rawString);
            emit newLog(qjsonError.errorString());
            continue;
        }

        processJSON(jdoc.object(), rawString);
    }
}

void PHD2::processJSON(const QJsonObject &jsonObj, QString rawString)
{
    PHD2MessageType messageType = PHD2_UNKNOWN;

    if (jsonObj.contains("Event"))
    {
        messageType = PHD2_EVENT;
        processPHD2Event(jsonObj);

        if (event == Alert)
            return;
    }
    else if (jsonObj.contains("error"))
    {
        messageType = PHD2_ERROR;
        processPHD2Error(jsonObj);
    }
    else if (jsonObj.contains("result"))
    {
        QJsonObject jsonResult = jsonObj["result"].toObject();
        if (jsonResult.contains("frame"))
        {
                messageType = PHD2_STAR_IMAGE;
                processStarImage(jsonResult);
        }
        else
        {
            messageType = PHD2_RESULT;
        }
    }

    if(messageType != PHD2_STAR_IMAGE)
        qCDebug(KSTARS_EKOS_GUIDE) << rawString;

    switch (connection)
    {
        case CONNECTING:
            if (event == Version)
                connection = CONNECTED;
            return;

        case CONNECTED:
            // If initial state is STOPPED, let us connect equipment
            if (state == STOPPED || state == PAUSED)
            {
                setEquipmentConnected(true);
            }
            else if (state == GUIDING || state == DITHERING)
            {
                connection = EQUIPMENT_CONNECTED;
                emit newStatus(Ekos::GUIDE_CONNECTED);
            }
            return;

        case DISCONNECTED:
            emit newStatus(Ekos::GUIDE_DISCONNECTED);
            break;

        case EQUIPMENT_CONNECTING:
            if (messageType == PHD2_RESULT)
            {
                connection = EQUIPMENT_CONNECTED;
                emit newStatus(Ekos::GUIDE_CONNECTED);
            }
            else if (messageType == PHD2_ERROR)
            {
                connection = EQUIPMENT_DISCONNECTED;
                emit newStatus(Ekos::GUIDE_DISCONNECTED);
            }
            return;

        case EQUIPMENT_CONNECTED:
        case EQUIPMENT_DISCONNECTED:
            break;

        case EQUIPMENT_DISCONNECTING:
            connection = EQUIPMENT_DISCONNECTED;
            //emit disconnected();
            return;
    }

    switch (state)
    {
        case GUIDING:
            break;

        case PAUSED:
            break;

        case STOPPED:
            break;

        default:
            break;
    }
}

void PHD2::processPHD2Event(const QJsonObject &jsonEvent)
{
    QString eventName = jsonEvent["Event"].toString();

    if (events.contains(eventName) == false)
    {
        emit newLog(i18n("Unknown PHD2 event: %1", eventName));
        return;
    }

    event = events.value(eventName);

    switch (event)
    {
        case Version:
            emit newLog(i18n("PHD2: Version %1", jsonEvent["PHDVersion"].toString()));
            break;

        case CalibrationComplete:
            //state = CALIBRATION_SUCCESSFUL;
            // It goes immediately to guiding until PHD implements a calibration-only method
            state = GUIDING;
            emit newLog(i18n("PHD2: Calibration Complete."));
            //emit guideReady();
            emit newStatus(Ekos::GUIDE_CALIBRATION_SUCESS);
            break;

        case StartGuiding:
            state = GUIDING;
            if (connection != EQUIPMENT_CONNECTED)
            {
                setConnectedRetries = 0;
                connection = EQUIPMENT_CONNECTED;
                emit newStatus(Ekos::GUIDE_CONNECTED);
            }
            emit newLog(i18n("PHD2: Guiding Started."));
            emit newStatus(Ekos::GUIDE_GUIDING);
            break;

        case Paused:
            state = PAUSED;
            emit newLog(i18n("PHD2: Guiding Paused."));
            emit newStatus(Ekos::GUIDE_SUSPENDED);
            break;

        case StartCalibration:
            state = CALIBRATING;
            emit newLog(i18n("PHD2: Calibration Started."));
            emit newStatus(Ekos::GUIDE_CALIBRATING);
            break;

        case AppState:
            processPHD2State(jsonEvent["State"].toString());
            break;

        case CalibrationFailed:
            state = CALIBRATION_FAILED;
            emit newLog(i18n("PHD2: Calibration Failed (%1).", jsonEvent["Reason"].toString()));
            emit newStatus(Ekos::GUIDE_CALIBRATION_ERROR);
            break;

        case CalibrationDataFlipped:
            emit newLog(i18n("Calibration Data Flipped."));
            break;

        case LoopingExposures:
            //emit newLog(i18n("PHD2: Looping Exposures."));
            break;

        case LoopingExposuresStopped:
            emit newLog(i18n("PHD2: Looping Exposures Stopped."));
            break;

        case Settling:        
        case SettleBegin:
        break;

        case SettleDone:
        {
            bool error = false;

            if (jsonEvent["Status"].toInt() != 0)
            {
                error = true;
                emit newLog(i18n("PHD2: Settling failed (%1).", jsonEvent["Error"].toString()));
            }

            if (state == GUIDING)
            {
                if (error)
                    state = STOPPED;
            }
            else if (state == DITHERING)
            {
                if (error)
                {
                    state = DITHER_FAILED;
                    //emit ditherFailed();
                    emit newStatus(GUIDE_DITHERING_ERROR);
                }
                else
                {
                    state = DITHER_SUCCESSFUL;
                    emit newStatus(Ekos::GUIDE_DITHERING_SUCCESS);
                }
            }
        }
        break;

        case StarSelected:
            emit newLog(i18n("PHD2: Star Selected."));
            break;

        case StarLost:
            emit newLog(i18n("PHD2: Star Lost."));
            emit newStatus(Ekos::GUIDE_ABORTED);
            break;

        case GuidingStopped:
            emit newLog(i18n("PHD2: Guiding Stopped."));
            state = STOPPED;
            //emit autoGuidingToggled(false);
            emit newStatus(Ekos::GUIDE_IDLE);
            break;

        case Resumed:
            emit newLog(i18n("PHD2: Guiding Resumed."));
            emit newStatus(Ekos::GUIDE_GUIDING);
            state = GUIDING;
            break;

        case GuideStep:
        {
            double diff_ra_pixels, diff_de_pixels, diff_ra_arcsecs, diff_de_arcsecs;
            diff_ra_pixels = jsonEvent["RADistanceRaw"].toDouble();
            diff_de_pixels = jsonEvent["DECDistanceRaw"].toDouble();

            diff_ra_arcsecs = 206.26480624709 * diff_ra_pixels * ccdPixelSizeX / mountFocalLength;
            diff_de_arcsecs = 206.26480624709 * diff_de_pixels * ccdPixelSizeY / mountFocalLength;

            emit newAxisDelta(diff_ra_arcsecs, diff_de_arcsecs);

            QJsonArray args2;
            args2<<32; //This is the size of the image that is being requested  32 x 32 pixels.
            sendJSONRPCRequest("get_star_image", args2); //This requests a star image for the guide view.
        }
        break;

        case GuidingDithered:
            emit newLog(i18n("PHD2: Guide Dithering."));
            state = DITHERING;
            emit newStatus(Ekos::GUIDE_DITHERING);
            break;

        case LockPositionSet:
            emit newLog(i18n("PHD2: Lock Position Set."));
            break;

        case LockPositionLost:
            emit newLog(i18n("PHD2: Lock Position Lost."));
            if (state == CALIBRATING)
                emit newStatus(Ekos::GUIDE_CALIBRATION_ERROR);
            break;

        case Alert:
            emit newLog(i18n("PHD2 %1: %2", jsonEvent["Type"].toString(), jsonEvent["Msg"].toString()));
            break;
    }
}

void PHD2::processStarImage(const QJsonObject &jsonStarFrame)
{
    //The width and height of the recieved PHD2 Star Image
   int width =  jsonStarFrame["width"].toInt();
   int height = jsonStarFrame["height"].toInt();

   //This sets up the Temp file which will be reused for subsequent captures
   QString filename = KSPaths::writableLocation(QStandardPaths::TempLocation) + QLatin1Literal("phd2.fits");

   //This section sets up the FITS File
   fitsfile *fptr;
   int status=0;
   long  fpixel = 1, naxis = 2, nelements, exposure;
   long naxes[2] = { width, height };
   fits_create_file(&fptr, QString("!"+filename).toLatin1().data(), &status);
   fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
    //Note, this is made up.  If you want the actual exposure time, you have to request it from PHD2
   exposure = 1;
   fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,"Total Exposure Time", &status);

   //This section takes the Pixels from the JSON Document
   //Then it converts from base64 to a QByteArray
   //Then it creates a datastream from the QByteArray to the pixel array for the FITS File   
   QByteArray converted = QByteArray::fromBase64(jsonStarFrame["pixels"].toString().toLocal8Bit());

   //This finishes up and closes the FITS file
   nelements = naxes[0] * naxes[1];
   fits_write_img(fptr, TUSHORT, fpixel, nelements, converted.data(), &status);
   fits_close_file(fptr, &status);
   fits_report_error(stderr, status);

   //This loads the FITS file in the Guide FITSView
   //Then it updates the Summary Screen
   bool imageLoad = guideFrame->loadFITS(filename, true);
   if (imageLoad)
   {
       guideFrame->updateFrame();
       guideFrame->setTrackingBox(QRect(0,0,width,height));
       emit newStarPixmap(guideFrame->getTrackingBoxPixmap());
   }

}

void PHD2::setGuideView(FITSView *guideView)
{
    guideFrame = guideView;
}

void PHD2::processPHD2State(const QString &phd2State)
{
    if (phd2State == "Stopped")
        state = STOPPED;
    else if (phd2State == "Selected")
        state = SELECTED;
    else if (phd2State == "Calibrating")
        state = CALIBRATING;
    else if (phd2State == "GUIDING")
        state = GUIDING;
    else if (phd2State == "LostLock")
        state = LOSTLOCK;
    else if (phd2State == "Paused")
        state = PAUSED;
    else if (phd2State == "Looping")
        state = LOOPING;
}

void PHD2::processPHD2Error(const QJsonObject &jsonError)
{
    QJsonObject jsonErrorObject = jsonError["error"].toObject();

    emit newLog(i18n("PHD2 Error: %1", jsonErrorObject["message"].toString()));

    if (state == DITHERING)
    {
        state = DITHER_FAILED;
        //emit ditherFailed();
        emit newStatus(GUIDE_DITHERING_ERROR);

        if (Options::ditherFailAbortsAutoGuide())
        {
            state = STOPPED;
            emit newStatus(GUIDE_ABORTED);
        }
        else
        {
            resume();
        }
     }
    /* switch (connection)
     {
         case CONNECTING:
         case CONNECTED:
             emit disconnected();
         break;

         default:
             break;
     }*/
}

void PHD2::sendJSONRPCRequest(const QString &method, const QJsonArray args)
{
    QJsonObject jsonRPC;

    jsonRPC.insert("jsonrpc", "2.0");
    jsonRPC.insert("method", method);

    if (args.empty() == false)
        jsonRPC.insert("params", args);

    jsonRPC.insert("id", methodID++);

    QJsonDocument json_doc(jsonRPC);

    //emit newLog(json_doc.toJson(QJsonDocument::Compact));
    qCDebug(KSTARS_EKOS_GUIDE) << json_doc.toJson(QJsonDocument::Compact);

    tcpSocket->write(json_doc.toJson(QJsonDocument::Compact));
    tcpSocket->write("\r\n");
}

void PHD2::setEquipmentConnected(bool enable)
{
    if (setConnectedRetries++ > MAX_SET_CONNECTED_RETRIES)
    {
        setConnectedRetries = 0;
        connection = EQUIPMENT_DISCONNECTED;
        emit newStatus(Ekos::GUIDE_DISCONNECTED);
        return;
    }

    if ((connection == EQUIPMENT_CONNECTED && enable == true) ||
        (connection == EQUIPMENT_DISCONNECTED && enable == false))
        return;

    if (enable)
        connection = EQUIPMENT_CONNECTING;
    else
        connection = EQUIPMENT_DISCONNECTING;

    QJsonArray args;

    // connected = enable
    args << enable;

    sendJSONRPCRequest("set_connected", args);
}

bool PHD2::calibrate()
{
    // We don't explicitly do calibration since it is done in the guide step by PHD2 anyway
    emit newStatus(Ekos::GUIDE_CALIBRATION_SUCESS);
    return true;
}

bool PHD2::guide()
{
    if (state == GUIDING)
    {
        emit newLog(i18n("PHD2: Guiding is already running."));
        emit newStatus(Ekos::GUIDE_GUIDING);
        return true;
    }

    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    QJsonArray args;
    QJsonObject settle;

    settle.insert("pixels", static_cast<double>(Options::ditherThreshold()));
    settle.insert("time", static_cast<int>(Options::ditherSettle()));
    settle.insert("timeout", static_cast<int>(Options::ditherTimeout()));

    // Settle param
    args << settle;
    // Recalibrate param
    args << false;

    sendJSONRPCRequest("guide", args);

    return true;
}

bool PHD2::abort()
{
    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    sendJSONRPCRequest("stop_capture");
    return true;
}

bool PHD2::suspend()
{
    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    QJsonArray args;

    // Paused param
    args << true;
    // FULL param
    args << "full";

    sendJSONRPCRequest("set_paused", args);

    return true;
}

bool PHD2::resume()
{
    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    QJsonArray args;

    // Paused param
    args << false;

    sendJSONRPCRequest("set_paused", args);

    return true;
}

bool PHD2::dither(double pixels)
{
    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    QJsonArray args;
    QJsonObject settle;

    settle.insert("pixels", static_cast<double>(Options::ditherThreshold()));
    settle.insert("time", static_cast<int>(Options::ditherSettle()));
    settle.insert("timeout", static_cast<int>(Options::ditherTimeout()));

    // Pixels
    args << pixels;
    // RA Only?
    args << false;
    // Settle
    args << settle;

    state = DITHERING;

    sendJSONRPCRequest("dither", args);

    return true;
}

bool PHD2::clearCalibration()
{
    if (connection != EQUIPMENT_CONNECTED)
    {
        emit newLog(i18n("PHD2 Error: Equipment not connected."));
        return false;
    }

    QJsonArray args;
    //This instructs PHD2 which calibration to clear.
    args << "mount";
    sendJSONRPCRequest("clear_calibration", args);

    return true;
}

}
