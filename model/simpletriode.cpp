#include "simpletriode.h"

struct SimpleTriodeResidual {
    SimpleTriodeResidual(double va, double vg, double ia) : va_(va), vg_(vg), ia_(ia) {}

    template <typename T>
    bool operator()(const T* const kg, const T* const vct, const T* const a, const T* const mu, T* residual) const {
        T e1t = va_ / mu[0] + vg_ + vct[0];
        if (e1t < 0.0) {
            e1t = mu[0] - mu[0];
        }
        T ia = pow(e1t, a[0]) / kg[0];
        residual[0] = ia_ - ia;
        return !(isnan(ia) || isinf(ia));
    }

private:
    const double va_;
    const double vg_;
    const double ia_;
};

SimpleTriode::SimpleTriode()
{
    parameter[TRI_KG] = new Parameter("Kg:", 0.7);
    parameter[TRI_VCT] = new Parameter("Vct:", 0.1);
    parameter[TRI_ALPHA] = new Parameter("Alpha:", 1.5);
    parameter[TRI_MU] = new Parameter("Mu:", 100.0);
}

void SimpleTriode::addSample(double va, double ia, double vg1, double vg2)
{
    problem.AddResidualBlock(
        new AutoDiffCostFunction<SimpleTriodeResidual, 1, 1, 1, 1, 1>(
            new SimpleTriodeResidual(va, vg1, ia)),
        NULL,
        parameter[TRI_KG]->getPointer(),
        parameter[TRI_VCT]->getPointer(),
        parameter[TRI_ALPHA]->getPointer(),
        parameter[TRI_MU]->getPointer());
}

double SimpleTriode::anodeCurrent(double va, double vg1, double vg2)
{
    double ia = 0.0;

    double e1t = va / parameter[TRI_MU]->getValue() + vg1 + parameter[TRI_VCT]->getValue();
    if (e1t > 0) {
        ia = pow(e1t, parameter[TRI_ALPHA]->getValue()) / parameter[TRI_KG]->getValue();
    }

    return ia;
}

void SimpleTriode::fromJson(QJsonObject source)
{
    if (source.contains("kg") && source["kg"].isDouble()) {
        parameter[TRI_KG]->setValue(source["kg"].toDouble());
    }

    if (source.contains("mu") && source["mu"].isDouble()) {
        parameter[TRI_MU]->setValue(source["mu"].toDouble());
    }

    if (source.contains("alpha") && source["alpha"].isDouble()) {
        parameter[TRI_ALPHA]->setValue(source["alpha"].toDouble());
    }

    if (source.contains("vct") && source["vct"].isDouble()) {
        parameter[TRI_VCT]->setValue(source["vct"].toDouble());
    }
}

void SimpleTriode::toJson(QJsonObject &destination, double vg1Max, double vg2Max)
{
    QJsonObject model;
    model["kg"] = parameter[TRI_KG]->getValue();
    model["mu"] = parameter[TRI_MU]->getValue();
    model["alpha"] = parameter[TRI_ALPHA]->getValue();
    model["vct"] = parameter[TRI_VCT]->getValue();

    QJsonObject triode;
    triode["vg1Max"] = vg1Max;
    triode["simple"] = model;

    destination["triode"] = triode;
}

void SimpleTriode::updateUI(QLabel *labels[], QLineEdit *values[])
{
    int i = 0;

    updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
    updateParameter(labels[i], values[i], parameter[TRI_VCT]); i++;
}

QString SimpleTriode::getName()
{
    return QString("Simple");
}

void SimpleTriode::setKg(double kg)
{
    parameter[TRI_KG]->setValue(kg);
}

void SimpleTriode::setMu(double mu)
{
    parameter[TRI_MU]->setValue(mu);
}

void SimpleTriode::setAlpha(double alpha)
{
    parameter[TRI_ALPHA]->setValue(alpha);
}

void SimpleTriode::setVct(double vct)
{
    parameter[TRI_VCT]->setValue(vct);
}

void SimpleTriode::setOptions()
{
    options.max_num_iterations = 100;
    options.minimizer_progress_to_stdout = true;

    setLowerBound(parameter[TRI_KG], 0.0000001); // Kg > 0
    setLimits(parameter[TRI_ALPHA], 1.0, 2.0); // 1.0 <= alpha <= 2.0
    setLimits(parameter[TRI_MU], 1.0, 1000.0); // 1.0 <= mu <= 1000.0
    setLimits(parameter[TRI_VCT], -2.0, 2.0); // -2.0 <= Vct <= 2.0
}
