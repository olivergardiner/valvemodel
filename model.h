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

/**
 * @brief The eTriodeParameter enum
 *
 * Defines indexes into the array of model Parameters. This is useful because the models considered
 * represent progressive refinements, i.e. the more complex models reuse parameters from the simpler
 * models and complex models only need to define the additional parameters that they add to the model
 */
enum eTriodeParameter {
    TRI_KG,
    TRI_KP,
    TRI_KVB,
    TRI_KVB2,
    TRI_VCT,
    TRI_ALPHA,
    TRI_MU
};

enum eModelType {
    SIMPLE_TRIODE,
    KOREN_TRIODE,
    IMPROVED_KOREN_TRIODE,
    KOREN_PENTODE,
    DERK_PENTODE,
    DERKE_PENTODE
};

/**
 * @brief sgn
 * @param val The value for which to compute the Signum
 * @return The Signum function of val
 *
 * Returns +1 if val > 0, -1 if val < 0 and 0 if val == 0
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief The Model class
 *
 * Model is the abstract base class for all of the device models and defines the public API as
 * well as providing some core methods that are common. The model uses the Ceres Solver library
 * to do a least squares fit of the underlying model.
 */
class Model : public UIBridge
{
public:
    /**
     * @brief addSample adds a sample to the sample set that will be use to fit the model
     * @param va The anode voltage
     * @param ia The anode current in mA
     * @param vg1 The grid voltage
     * @param vg2 For pentodes only, the screen grid voltage
     */
	virtual void addSample(double va, double ia, double vg1, double vg2 = 0.0) = 0;
    /**
     * @brief fromJson reads the model parameters from a Json object
     * @param source The Json object to read
     */
    virtual void fromJson(QJsonObject source) = 0;
    /**
     * @brief toJson writes the current model parameters to a Json object
     * @param destination The Json object to write to
     */
    virtual void toJson(QJsonObject destination) = 0;
    /**
     * @brief anodeCurrent calculates the modelled anode current
     * @param va The anode voltage
     * @param vg1 The grid voltage
     * @param vg2 For pentodes only, the screen grid voltage
     * @return The anode current in mA
     */
    virtual double anodeCurrent(double va, double vg1, double vg2 = 0.0) = 0;
    virtual double anodeVoltage(double ia, double vg1, double vg2 = 0.0);
    virtual QString getName() = 0;

    void solve();

 protected:
    /**
     * @brief problem The Ceres Problem used for model fitting
     */
	Problem problem;
    /**
     * @brief parameter The array of 8 model Parameters linked to the UI
     */
	Parameter *parameter[8];
    /**
     * @brief options The options to be used by Ceres for solving the model approximation
     */
    Solver::Options options;

    void setLowerBound(Parameter* parameter, double lowerBound);
    void setUpperBound(Parameter* parameter, double upperBound);
    void setLimits(Parameter* parameter, double lowerBound, double upperBound);
    virtual void setOptions() = 0;
    double korenCurrent(double va, double vg, double kp, double kvb, double a, double mu);
    double improvedKorenCurrent(double va, double vg, double kp, double kvb, double kvb2, double vct, double a, double mu);
};

