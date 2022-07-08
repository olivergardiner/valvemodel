#include "korentriode.h"

struct KorenTriodeResidual {
    KorenTriodeResidual(double va, double vg, double ia) : va_(va), vg_(vg), ia_(ia) {}

    template <typename T>
    bool operator()(const T* const kg, const T* const kp, const T* const kvb, const T* const a, const T* const mu, T* residual) const {
        T e1t = log(1.0 + exp(kp[0] * (1.0 / mu[0] + vg_ / sqrt(kvb[0] + va_ * va_))));
        T ia = pow((va_ / kp[0]) * e1t, a[0]) / kg[0];
        residual[0] = ia_ - ia;
        return !(isnan(ia) || isinf(ia));
    }

private:
    const double va_;
    const double vg_;
    const double ia_;
};

KorenTriode::KorenTriode()
{
    parameter[TRI_KP] = new Parameter("Kp:", 500.0);
    parameter[TRI_KVB] = new Parameter("Kvb:", 300.0);
}

void KorenTriode::addSample(double va, double ia, double vg1, double vg2)
{
    problem.AddResidualBlock(
        new AutoDiffCostFunction<KorenTriodeResidual, 1, 1, 1, 1, 1, 1>(
            new KorenTriodeResidual(va, vg1, ia)),
        NULL,
        parameter[TRI_KG]->getPointer(),
        parameter[TRI_KP]->getPointer(),
        parameter[TRI_KVB]->getPointer(),
        parameter[TRI_ALPHA]->getPointer(),
        parameter[TRI_MU]->getPointer());
}

double KorenTriode::anodeCurrent(double va, double vg1, double vg2)
{
    return korenCurrent(va, vg1,
        parameter[TRI_KP]->getValue(),
        parameter[TRI_KVB]->getValue(),
        parameter[TRI_ALPHA]->getValue(),
        parameter[TRI_MU]->getValue()) / parameter[TRI_KG]->getValue();

}

void KorenTriode::fromJson(QJsonObject source)
{
    SimpleTriode::fromJson(source);

    if (source.contains("kp") && source["kp"].isDouble()) {
        parameter[TRI_KP]->setValue(source["kp"].toDouble());
    }

    if (source.contains("kvb") && source["kvb"].isDouble()) {
        parameter[TRI_KVB]->setValue(source["kvb"].toDouble());
    }
}

void KorenTriode::toJson(QJsonObject &destination, double vg1Max, double vg2Max)
{
    QJsonObject model;
    model["kg"] = parameter[TRI_KG]->getValue();
    model["mu"] = parameter[TRI_MU]->getValue();
    model["alpha"] = parameter[TRI_ALPHA]->getValue();
    model["vct"] = parameter[TRI_VCT]->getValue();
    model["kp"] = parameter[TRI_KP]->getValue();
    model["kvb"] = parameter[TRI_KVB]->getValue();

    QJsonObject triode;
    triode["vg1Max"] = vg1Max;
    triode["koren"] = model;

    destination["triode"] = triode;
}

void KorenTriode::updateUI(QLabel *labels[], QLineEdit *values[])
{
    int i = 0;

    updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KP]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KVB]); i++;
}

QString KorenTriode::getName()
{
    return QString("Koren");
}

void KorenTriode::setOptions()
{
    SimpleTriode::setOptions();

    setLimits(parameter[TRI_KVB], 0.0, 10000.0); // 0 <= Kvb <= 10000.0
    options.linear_solver_type = ceres::CGNR;
    options.preconditioner_type = ceres::JACOBI;
}
