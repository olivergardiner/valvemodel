#ifndef UIBRIDGE_H
#define UIBRIDGE_H

#include <QLineEdit>
#include <QLabel>

#include "parameter.h"

class UIBridge
{
public:
    UIBridge();

    virtual void updateUI(QLabel *labels[], QLineEdit *values[]) = 0;

protected:
    void updateParameter(QLabel *uiLabel, QLineEdit *uiValue, QString label, double value);
    void updateParameter(QLabel *uiLabel, QLineEdit *uiValue, Parameter *parameter);
};

#endif // UIBRIDGE_H
