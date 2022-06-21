#include "modelfactory.h"

Model *ModelFactory::createModel(int modelType)
{
    switch ((eModelType) modelType) {
    case SIMPLE_TRIODE:
        return new SimpleTriode();
    case KOREN_TRIODE:
        return new KorenTriode();
    case IMPROVED_KOREN_TRIODE:
        return new ImprovedKorenTriode();
    case KOREN_PENTODE:
        break;
    case DERK_PENTODE:
        break;
    case DERKE_PENTODE:
        break;
    }

    return nullptr;
}
