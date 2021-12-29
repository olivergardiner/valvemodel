#pragma once

#include "korentriode.h"

class ImprovedKorenTriode : public KorenTriode
{
public:
    ImprovedKorenTriode();

	virtual void addSample(double va, double ia, double vg1, double vg2 = 0.0);
	virtual double anodeCurrent(double va, double vg1, double vg2 = 0.0);
    virtual void fromJson(QJsonObject source);
	virtual void toJson(QJsonObject destination);
    virtual void updateUI(QLabel *labels[], QLineEdit *values[]);
    virtual QString getName();

protected:
	void setOptions();
};
