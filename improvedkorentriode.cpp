#include "improvedkorentriode.h"

struct ImprovedKorenTriodeResidual {
    ImprovedKorenTriodeResidual(double va, double vg, double ia) : va_(va), vg_(vg), ia_(ia) {}

    template <typename T>
    bool operator()(const T* const kg, const T* const kp, const T* const kvb, const T* const kvb2, const T* const vct, const T* const a, const T* const mu, T* residual) const {
        T e2t = log(1.0 + exp(kp[0] * (1.0 / mu[0] + (vg_ + vct[0]) / sqrt(kvb[0] + va_ * va_ + kvb2[0] * va_))));
        T ia = pow((va_ / kp[0]) * e2t, a[0]) / kg[0];
        residual[0] = ia_ - ia;
        return !(isnan(ia) || isinf(ia));
    }

private:
    const double va_;
    const double vg_;
    const double ia_;
};

ImprovedKorenTriode::ImprovedKorenTriode()
{
    parameter[TRI_KVB2] = new Parameter("Kvb2:", 30.0);
}

void ImprovedKorenTriode::addSample(double va, double ia, double vg1, double vg2)
{
    problem.AddResidualBlock(
        new AutoDiffCostFunction<ImprovedKorenTriodeResidual, 1, 1, 1, 1, 1, 1, 1, 1>(
            new ImprovedKorenTriodeResidual(va, vg1, ia)),
        NULL,
        parameter[TRI_KG]->getPointer(),
        parameter[TRI_KP]->getPointer(),
        parameter[TRI_KVB]->getPointer(),
        parameter[TRI_KVB2]->getPointer(),
        parameter[TRI_VCT]->getPointer(),
        parameter[TRI_ALPHA]->getPointer(),
        parameter[TRI_MU]->getPointer());
}

double ImprovedKorenTriode::anodeCurrent(double va, double vg1, double vg2)
{
    return improvedKorenCurrent(va, vg1,
        parameter[TRI_KP]->getValue(),
        parameter[TRI_KVB]->getValue(),
        parameter[TRI_KVB2]->getValue(),
        parameter[TRI_VCT]->getValue(),
        parameter[TRI_ALPHA]->getValue(),
        parameter[TRI_MU]->getValue()) / parameter[TRI_KG]->getValue();
}

void ImprovedKorenTriode::fromJson(QJsonObject source)
{
    KorenTriode::fromJson(source);

    if (source.contains("kvb2") && source["kvb2"].isDouble()) {
        parameter[TRI_KVB2]->setValue(source["kvb2"].toDouble());
    }
}

void ImprovedKorenTriode::toJson(QJsonObject destination)
{
}

void ImprovedKorenTriode::updateUI(QLabel *labels[], QLineEdit *values[])
{
    int i = 0;

    updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KP]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KVB]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KVB2]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_VCT]); i++;
}

QString ImprovedKorenTriode::getName()
{
    return QString("Improved Koren");
}

void ImprovedKorenTriode::setOptions()
{
    KorenTriode::setOptions();

    setLimits(parameter[TRI_KVB2], 0.0, 1000.0); // 0.0 <= Kvb2 <= 1000.0
    setLimits(parameter[TRI_VCT], 0.0, 2.0); // 0.0 <= Vct <= 2.0
    options.linear_solver_type = ceres::CGNR;
    options.preconditioner_type = ceres::JACOBI;
}
