
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
    function<const double (const Path&)> control_variate;
    double control_variate_mean = 0.0;
    bool m_PathwiseDifferentiation = false;
public:
    MCModel(Option &option, double S0, double sigma, double r, double h = 0.01, SDESolver sdeSolver = Explicit): Model(option, S0, sigma, r) {
        this->m_h = h;
        this->m_Solver = sdeSolver;
    }

    void define_control_variate(function<const double (const Path&)> control_variate, double control_variate_mean) {
        this->control_variate = control_variate;
        this->control_variate_mean = control_variate_mean;
    }

    ModelResult simulate(Simulator simulator, int simulations, int path_size = 1, bool pathwise = false) {
        clock_t start = clock();
        this->simulation_vector.clear();
        this->simulation_vector = simulator.simulate(*this, simulations, path_size);
        this->m_PathwiseDifferentiation = pathwise;
        auto price_and_variance = this->calcPrice();
        double delta = this->calcDelta();
        double gamma = this->calcGamma();
        double vega  = this->calcVega();
        ModelResult result;
        result.setPrice(price_and_variance.first);
        result.setPriceVariance(price_and_variance.second);
        result.setDelta(delta);
        result.setGamma(gamma);
        result.setVega(vega);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
    pair<double,double> calcPrice() const override;
    double calcDelta() const override;
    double calcGamma() const override;
    double calcVega() const override;
};


#endif //SIMULATIONMETHODS_MCMODEL_H
