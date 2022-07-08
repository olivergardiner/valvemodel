#pragma once

#include "simpletriode.h"
#include "korentriode.h"
#include "improvedkorentriode.h"

class ModelFactory
{
public:
    static Model *createModel(int modelType);
};

