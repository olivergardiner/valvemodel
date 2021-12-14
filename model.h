#pragma once

#include <QJsonObject>
#include <QString>

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "parameter.h"
#include "uibridge.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

enum eTriodeParameter {
    TRI_KG,
    TRI_KP,
    TRI_KVB,
    TRI_KVB2,
    TRI_VCT,
    TRI_ALPHA,
    TRI_MU
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class Model : public UIBridge
{
public:
	virtual void addSample(double va, double ia, double vg1, double vg2 = 0.0) = 0;
    virtual void fromJson(QJsonObject source) = 0;
    virtual void toJson(QJsonObject destination) = 0;
    virtual double anodeCurrent(double va, double vg1, double vg2 = 0.0) = 0;
    virtual double anodeVoltage(double ia, double vg1, double vg2 = 0.0);
    virtual QString getName() = 0;

    void solve();

 protected:
	Problem problem;
	Parameter *parameter[8];
    Solver::Options options;

    void setLowerBound(Parameter* parameter, double lowerBound);
    void setUpperBound(Parameter* parameter, double upperBound);
    void setLimits(Parameter* parameter, double lowerBound, double upperBound);
    virtual void setOptions() = 0;
    double korenCurrent(double va, double vg, double kp, double kvb, double a, double mu);
    double improvedKorenCurrent(double va, double vg, double kp, double kvb, double kvb2, double vct, double a, double mu);
};

