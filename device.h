#pragma once

#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QJsonDocument>
#include <QJsonValue>
#include <QJsonArray>
#include <QJsonObject>
#include <QGraphicsItemGroup>

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "parameter.h"
#include "model.h"
#include "uibridge.h"
#include "plot.h"
#include "simpletriode.h"
#include "korentriode.h"
#include "improvedkorentriode.h"

enum eModelType {
    SIMPLE_TRIODE,
    KOREN_TRIODE,
    IMPROVED_KOREN_TRIODE,
    KOREN_PENTODE,
    DERK_PENTODE,
    DERKE_PENTODE
};

enum eModelDeviceType {
    MODEL_TRIODE,
    MODEL_PENTODE
};

enum ePlotType {
    PLOT_TRIODE_ANODE,
    PLOT_TRIODE_TRANSFER,
    PLOT_PENTODE_ANODE,
    PLOT_PENTODE_TRANSFER,
    PLOT_PENTODE_SCREEN
};

class Device : UIBridge
{
public:
    Device(int _modelDeviceType);
    Device(QJsonDocument model);

    double getParameter(int index) const;

    void solve();

    double anodeCurrent(double va, double vg1, double vg2 = 0);
    double anodeVoltage(double ia, double vg1, double vg2 = 0);

    void updateUI(QLabel *labels[], QLineEdit *values[]);
    void updateModelSelect(QComboBox *select);
    void selectModel(int index);
    void anodeAxes(Plot *plot);
    void transferAxes(Plot *plot);
    QGraphicsItemGroup *anodePlot(Plot *plot);
    QGraphicsItemGroup *transferPlot(Plot *plot);
    double interval(double maxValue);

    int getModelType() const;
    void setModelType(int newModelType);

    double getVaMax() const;
    double getIaMax() const;
    double getVg1Max() const;
    double getVg2Max() const;
    double getPaMax() const;

    int getDeviceType() const;

    void setDeviceType(int newDeviceType);

    QString getName();

private:
    int deviceType = MODEL_TRIODE;
    int modelType = IMPROVED_KOREN_TRIODE;

    QList<Model *> models;
    Model *currentModel = nullptr;

    QString name;

    double vaMax;
    double iaMax;
    double vg1Max;
    double vg2Max;
    double ig2Max;
    double paMax;
};
