#include "model.h"

double Model::anodeVoltage(double ia, double vg1, double vg2)
{
    double va = 100.0;

    double iaTest = anodeCurrent(va, vg1, vg2);
    double gradient = iaTest - anodeCurrent(va - 1.0, vg1, vg2);
    double iaErr = ia - iaTest;

    while ((abs(iaErr) / ia) > 0.01) {
        va = va + iaErr / gradient;
        iaTest = anodeCurrent(va, vg1, vg2);
        gradient = iaTest - anodeCurrent(va - 1.0, vg1, vg2);
        iaErr = ia - iaTest;
    }

    return va;
}

void Model::solve()
{
    setOptions();

    Solver::Summary summary;
    Solve(options, &problem, &summary);

    qInfo(summary.BriefReport().c_str());
}

void Model::setLowerBound(Parameter* parameter, double lowerBound)
{
    problem.SetParameterLowerBound(parameter->getPointer(), 0, lowerBound);
}

void Model::setUpperBound(Parameter* parameter, double upperBound)
{
    problem.SetParameterUpperBound(parameter->getPointer(), 0, upperBound);
}

void Model::setLimits(Parameter* parameter, double lowerBound, double upperBound)
{
    problem.SetParameterLowerBound(parameter->getPointer(), 0, lowerBound);
    problem.SetParameterUpperBound(parameter->getPointer(), 0, upperBound);
}

double Model::korenCurrent(double va, double vg, double kp, double kvb, double a, double mu)
{
    double x1 = std::sqrt(kvb + va * va);
    double x2 = kp * (1 / mu + vg / x1);
    double x3 = std::log(1.0 + std::exp(x2));
    double et = (va / kp) * x3;

    if (et < 0.0) {
        et = 0.0;
    }

    return pow(et, a);
}

double Model::improvedKorenCurrent(double va, double vg, double kp, double kvb, double kvb2, double vct, double a, double mu)
{
    double x1 = std::sqrt(kvb + va * va + va * kvb2);
    double x2 = kp * (1 / mu + (vg + vct) / x1);
    double x3 = std::log(1.0 + std::exp(x2));
    double et = (va / kp) * x3;

    if (et < 0.0) {
        et = 0.0;
    }

    return pow(et, a);
}
