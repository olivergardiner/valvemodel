#ifndef DEVICEMODEL_H
#define DEVICEMODEL_H

#include <QLabel>
#include <QLineEdit>
#include <QJsonDocument>
#include <QJsonValue>
#include <QJsonArray>
#include <QJsonObject>

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

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

enum eTriodeParameter {
    TRI_KG,
    TRI_KP,
    TRI_KVB,
    TRI_KVB2,
    TRI_VCT,
    TRI_ALPHA,
    TRI_MU
};

class DeviceModel
{
public:
    DeviceModel(int _modelDeviceType);
    DeviceModel(QJsonDocument model);

    double getParameter(int index) const;

    void addTriodeSample(double va, double vg1, double ia);
    void addPentodeSample(double va, double vg1, double vg2, double ia);
    void solve();

    double anodeCurrent(double va, double vg1);
    double anodeCurrent(double va, double vg1, double vg2);

    void updateUI(QLabel *labels[], QLineEdit *values[]);

    int getModelType() const;
    void setModelType(int newModelType);

private:
    Problem problem;
    int modelDeviceType = MODEL_TRIODE;
    int modelType = IMPROVED_KOREN_TRIODE;

    // Model variables
    QString parameterName[8];
    double parameterValue[8];

    double korenCurrent(double va, double vg, double kp, double kvb, double a, double mu);
    double improvedKorenCurrent(double va, double vg, double kp, double kvb, double kvb2, double vct, double a, double mu);

    void updateParameter(QLabel *label, QLineEdit *value, QString name, double x);
};

#endif // DEVICEMODEL_H
