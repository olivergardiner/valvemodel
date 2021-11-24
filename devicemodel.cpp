#include "devicemodel.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
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

struct ImprovedKorenTriodeResidual {
    ImprovedKorenTriodeResidual(double va, double vg, double ia) : va_(va), vg_(vg), ia_(ia) {}

    template <typename T>
    bool operator()(const T* const kg, const T* const kp, const T* const kvb, const T* const kvb2, const T* const vct, const T* const a, const T* const mu, T* residual) const {
        T e2t = log(1.0 + exp(kp[0] * (1.0 / mu[0] + (vg_ + vct[0])/ sqrt(kvb[0] + va_ * va_ + kvb2[0] * va_))));
        T ia = pow((va_ / kp[0]) * e2t, a[0]) / kg[0];
        residual[0] = ia_ - ia;
        return !(isnan(ia) || isinf(ia));
    }

private:
    const double va_;
    const double vg_;
    const double ia_;
};

DeviceModel::DeviceModel(int _modelDeviceType) : modelDeviceType(_modelDeviceType)
{
    if (modelDeviceType == MODEL_TRIODE) {
        initTriode();
        modelType = IMPROVED_KOREN_TRIODE;
    }
}

DeviceModel::DeviceModel(QJsonDocument modelDocument)
{
    if (modelDocument.isObject()) {
        QJsonObject model = modelDocument.object();

        if (model.contains("triode") && model["triode"].isObject()) {
            QJsonObject triode = model["triode"].toObject();

            modelDeviceType = MODEL_TRIODE;
            initTriode();

            if (triode.contains("kg") && triode["kg"].isDouble()) {
                parameter[TRI_KG]->setValue(triode["kg"].toDouble());
            }

            if (triode.contains("mu") && triode["mu"].isDouble()) {
                parameter[TRI_MU]->setValue(triode["mu"].toDouble());
            }

            if (triode.contains("alpha") && triode["alpha"].isDouble()) {
                parameter[TRI_ALPHA]->setValue(triode["alpha"].toDouble());
            }

            if (triode.contains("kp") && triode["kp"].isDouble()) {
                parameter[TRI_KP]->setValue(triode["kp"].toDouble());
            }

            if (triode.contains("kvb") && triode["kvb"].isDouble()) {
                parameter[TRI_KVB]->setValue(triode["kvb"].toDouble());
            }

            if (triode.contains("kvb2") && triode["kvb2"].isDouble()) {
                parameter[TRI_KVB2]->setValue(triode["kvb2"].toDouble());
            }

            if (triode.contains("vct") && triode["vct"].isDouble()) {
                parameter[TRI_VCT]->setValue(triode["vct"].toDouble());
            }
        }
    }
}

void DeviceModel::initTriode()
{
    parameter[TRI_KG] = new Parameter("Kg:", 0.7);
    parameter[TRI_KP] = new Parameter("Kp:", 0.7);
    parameter[TRI_KVB] = new Parameter("Kvb:", 0.7);
    parameter[TRI_KVB2] = new Parameter("Kvb2:", 0.7);
    parameter[TRI_VCT] = new Parameter("Vct:", 0.7);
    parameter[TRI_ALPHA] = new Parameter("Alpha:", 0.7);
    parameter[TRI_MU] = new Parameter("Mu:", 100.0);
}

void DeviceModel::initPentode()
{

}

void DeviceModel::addTriodeSample(double va, double vg1, double ia)
{
    if (modelType == SIMPLE_TRIODE) {
        problem.AddResidualBlock(
            new AutoDiffCostFunction<SimpleTriodeResidual, 1, 1, 1, 1, 1>(
                new SimpleTriodeResidual(va, vg1, ia)),
            NULL,
            parameter[TRI_KG]->getPointer(),
            parameter[TRI_VCT]->getPointer(),
            parameter[TRI_ALPHA]->getPointer(),
            parameter[TRI_MU]->getPointer());
    } else if (modelType == KOREN_TRIODE) {
        problem.AddResidualBlock(
            new AutoDiffCostFunction<KorenTriodeResidual, 1, 1, 1, 1, 1, 1>(
                new KorenTriodeResidual(va, vg1, ia)),
            NULL,
            parameter[TRI_KG]->getPointer(),
            parameter[TRI_KP]->getPointer(),
            parameter[TRI_KVB]->getPointer(),
            parameter[TRI_ALPHA]->getPointer(),
            parameter[TRI_MU]->getPointer());
    } else if (modelType == IMPROVED_KOREN_TRIODE) {
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
}

void DeviceModel::addPentodeSample(double va, double vg1, double vg2, double ia)
{

}

void DeviceModel::solve()
{
    Solver::Options options;

    if (modelDeviceType == MODEL_TRIODE) {
        setLowerBound(parameter[TRI_KG], 0.0000001); // Kg > 0
        setLimits(parameter[TRI_ALPHA], 1.0, 2.0); // 1.0 <= alpha <= 2.0
        setLimits(parameter[TRI_MU], 1.0, 1000.0); // 1.0 <= mu <= 1000.0
        setLimits(parameter[TRI_VCT], -2.0, 2.0); // -2.0 <= Vct <= 2.0

        if (modelType == KOREN_TRIODE || modelType == IMPROVED_KOREN_TRIODE) {
            setLimits(parameter[TRI_KVB], 0.0, 10000.0); // 0 <= Kvb <= 10000.0
            options.linear_solver_type = ceres::CGNR;
            options.preconditioner_type = ceres::JACOBI;
        }

        if (modelType == IMPROVED_KOREN_TRIODE) {
            setLimits(parameter[TRI_KVB2], 0.0, 1000.0); // 0.0 <= Kvb2 <= 1000.0
            setLimits(parameter[TRI_VCT], 0.0, 2.0); // 0.0 <= Vct <= 2.0
            options.linear_solver_type = ceres::CGNR;
            options.preconditioner_type = ceres::JACOBI;
        }
    }

    options.max_num_iterations = 100;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    qInfo(summary.BriefReport().c_str());
}

void DeviceModel::setLowerBound(Parameter *parameter, double lowerBound)
{
    problem.SetParameterLowerBound(parameter->getPointer(), 0, lowerBound);
}

void DeviceModel::setUpperBound(Parameter *parameter, double upperBound)
{
    problem.SetParameterUpperBound(parameter->getPointer(), 0, upperBound);
}

void DeviceModel::setLimits(Parameter *parameter, double lowerBound, double upperBound)
{
    problem.SetParameterLowerBound(parameter->getPointer(), 0, lowerBound);
    problem.SetParameterUpperBound(parameter->getPointer(), 0, upperBound);
}

double DeviceModel::anodeCurrent(double va, double vg1)
{
    // Only va and vg1 supplied so necessarily a triode
    double ia = 0.0;

    if (modelType == SIMPLE_TRIODE) {
        double e1t = va / parameter[TRI_MU]->getValue() + vg1 + parameter[TRI_VCT]->getValue();
        if (e1t > 0) {
            ia = pow(e1t, parameter[TRI_ALPHA]->getValue()) / parameter[TRI_KG]->getValue();
        }
    } else if (modelType == KOREN_TRIODE) {
        ia = korenCurrent(va, vg1, parameter[TRI_KP]->getValue(), parameter[TRI_KVB]->getValue(), parameter[TRI_ALPHA]->getValue(), parameter[TRI_MU]->getValue()) / parameter[TRI_KG]->getValue();
    } else if (modelType == IMPROVED_KOREN_TRIODE) {
        ia = improvedKorenCurrent(va, vg1, parameter[TRI_KP]->getValue(), parameter[TRI_KVB]->getValue(), parameter[TRI_KVB2]->getValue(), parameter[TRI_VCT]->getValue(), parameter[TRI_ALPHA]->getValue(), parameter[TRI_MU]->getValue()) / parameter[TRI_KG]->getValue();
    }

    return ia;
}

double DeviceModel::anodeCurrent(double va, double vg1, double vg2)
{
    // Vg2 also provided so the pentode models should be being used
    double ia = 0.0;

    return ia;
}

double DeviceModel::korenCurrent(double va, double vg, double kp, double kvb, double a, double mu)
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

double DeviceModel::improvedKorenCurrent(double va, double vg, double kp, double kvb, double kvb2, double vct, double a, double mu)
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

void DeviceModel::updateUI(QLabel *labels[], QLineEdit *values[])
{
    int i = 0;

    if (modelType == SIMPLE_TRIODE) {
        updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_VCT]); i++;
    } else if (modelType == KOREN_TRIODE) {
        updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KP]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KVB]); i++;
    } else if (modelType == IMPROVED_KOREN_TRIODE) {
        updateParameter(labels[i], values[i], parameter[TRI_MU]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KG]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_ALPHA]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KP]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KVB]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_KVB2]); i++;
        updateParameter(labels[i], values[i], parameter[TRI_VCT]); i++;
    }
}

void DeviceModel::updateParameter(QLabel *label, QLineEdit *value, Parameter *parameter)
{
    label->setText(parameter->getName());
    label->setVisible(true);

    char number[32];

    sprintf(number, "%.3f", parameter->getValue());

    int length = strlen(number);
    for (int i=length-1;i >= 0; i--) {
        char test = number[i];
        if (test == '0' || test == '.') {
            number[i] = 0;
        }

        if (test != '0') {
            break;
        }
    }

    value->setText(number);
    value->setVisible(true);
}

int DeviceModel::getModelType() const
{
    return modelType;
}

void DeviceModel::setModelType(int newModelType)
{
    modelType = newModelType;
}

double DeviceModel::getParameter(int index) const
{
    return parameter[index]->getValue();
}

double DeviceModel::getVaMax() const
{
    return vaMax;
}

double DeviceModel::getIaMax() const
{
    return iaMax;
}

double DeviceModel::getVg1Max() const
{
    return vg1Max;
}

double DeviceModel::getVg2Max() const
{
    return vg2Max;
}

double DeviceModel::getPaMax() const
{
    return paMax;
}
