#ifndef SIMPLETRIODE_H
#define SIMPLETRIODE_H

#include "model.h"

class SimpleTriode : public Model
{
public:
    SimpleTriode();

	virtual void addSample(double va, double ia, double vg1, double vg2 = 0.0);
	virtual double anodeCurrent(double va, double vg1, double vg2 = 0.0);
    virtual void fromJson(QJsonObject source);
	virtual void toJson(QJsonObject destination);
    virtual void updateUI(QLabel *labels[], QLineEdit *values[]);
    virtual QString getName();

	void setKg(double kg);
	void setMu(double kg);
	void setAlpha(double kg);
	void setVct(double kg);

protected:
	void setOptions();
};

#endif // SIMPLETRIODE_H
