
#ifndef SIMULATIONMETHODS_MCMODEL_H
#define SIMULATIONMETHODS_MCMODEL_H

#include "Model.h"
#include "Simulator.h"
#include "ModelResult.h"
#include <vector>
#include <ctime>

using namespace std;

class MCModel: public Model {
protected:
    vector<Path> simulation_vector;
public:
    MCModel(Option option, double S0, double sigma, double r, double h = 0.01): Model(option, S0, sigma, r) {
        this->m_h = h;
    }
    ModelResult simulate(Simulator simulator, int simulations, int path_size = 1) {
        clock_t start = clock();
        this->simulation_vector = simulator.simulate(*this, simulations, path_size);
        double price = this->calcPrice();
        double delta = this->calcDelta();
        double gamma = this->calcGamma();
        double vega  = this->calcVega();
        ModelResult result;
        result.setPrice(price);
        result.setDelta(delta);
        result.setGamma(gamma);
        result.setVega(vega);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
    double calcPrice() const override;
    double calcDelta() const override;
    double calcGamma() const override;
    double calcVega() const override;
};


#endif //SIMULATIONMETHODS_MCMODEL_H
