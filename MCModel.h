
#ifndef SIMULATIONMETHODS_MCMODEL_H
#define SIMULATIONMETHODS_MCMODEL_H

#include "Model.h"

class MCModel: public Model {
public:
    MCModel(Option option, double S0, double sigma, double r): Model(option, S0, sigma, r) {}
};


#endif //SIMULATIONMETHODS_MCMODEL_H
