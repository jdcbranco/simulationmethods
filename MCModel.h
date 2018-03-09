
#ifndef SIMULATIONMETHODS_MCMODEL_H
#define SIMULATIONMETHODS_MCMODEL_H

#include "ModelResult.h"
#include "Model.h"
#include "Simulator.h"
#include "SensitivityModel.h"

#include <vector>
#include <ctime>
#include <algorithm>
#include <numeric>

using namespace std;

template<class OPTION>
class MCModel: public Model<OPTION> {
    using Model<OPTION>::m_Option;
    using Model<OPTION>::m_S0;
    using Model<OPTION>::m_Sigma;
    using Model<OPTION>::discount;
    using Model<OPTION>::m_h;
protected:
    vector<Path> simulation_vector;
    function<const double (const Path&)> control_variate;
    double control_variate_mean = 0.0;
    SensitivityMethod m_SensitivityMethod = SensitivityMethod::FiniteDifference;
public:
    MCModel(OPTION &option, double S0, double sigma, double r, double h = 0.01, SDESolver sdeSolver = Explicit): Model<OPTION>(option, S0, sigma, r) {
        this->m_h = h;
        this->m_Solver = sdeSolver;
    }

    void define_control_variate(function<const double (const Path&)> control_variate, double control_variate_mean) {
        this->control_variate = control_variate;
        this->control_variate_mean = control_variate_mean;
    }

    ModelResult simulate(Simulator simulator, SensitivityModel<OPTION> &sensitivityModel, int simulations, int path_size = 1) {
        clock_t start = clock();
        this->simulation_vector.clear();
        this->simulation_vector = simulator.simulate(*this, simulations, path_size);
        auto price = this->calcPrice();
        auto delta = this->calcDelta(sensitivityModel);
        auto gamma = this->calcGamma(sensitivityModel);
        auto vega  = this->calcVega(sensitivityModel);
        ModelResult result;
        result.setModelType(ModelType::MonteCarlo);
        result.setDeltaMethod(delta.second);
        result.setGammaMethod(gamma.second);
        result.setVegaMethod(vega.second);
        result.setPrice(price.first);
        result.setPriceVariance(price.second);
        result.setDelta(delta.first);
        result.setGamma(gamma.first);
        result.setVega(vega.first);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }

    pair<double,double> calcPrice() const override;
    pair<double,SensitivityMethod> calcDelta() const override {
        Greeks_by_FD::CentralDifferencesSensitivityModel<OPTION> fd_method(*this,this->m_h);
        return calcDelta(fd_method);
    };
    pair<double,SensitivityMethod> calcGamma() const override {
        Greeks_by_FD::CentralDifferencesSensitivityModel<OPTION> fd_method(*this,this->m_h);
        return calcGamma(fd_method);
    };
    pair<double,SensitivityMethod> calcVega() const override {
        Greeks_by_FD::CentralDifferencesSensitivityModel<OPTION> fd_method(*this,this->m_h);
        return calcVega(fd_method);
    };
    pair<double,SensitivityMethod> calcDelta(SensitivityModel<OPTION> &sensitivityModel) const;
    pair<double,SensitivityMethod> calcGamma(SensitivityModel<OPTION> &sensitivityModel) const;
    pair<double,SensitivityMethod> calcVega(SensitivityModel<OPTION> &sensitivityModel) const;
};

template<class OPTION> pair<double,double> MCModel<OPTION>::calcPrice() const {
    if(!control_variate) {
        double sum = 0.0, sum_of_squares = 0.0;
        for(auto path: simulation_vector) {
            double discounted = discount(m_Option.payoff(path));
            sum += discounted;
            sum_of_squares += discounted * discounted;
        }
        double size = simulation_vector.size();
        double estimate = size > 0 ? sum / size : NAN;
        double variance = ((sum_of_squares / size) - estimate * estimate) / size;
        return pair<double, double>(estimate, variance);
    } else {
        cout << "* Control variate version *" << endl;
        vector<double> payoffs, control_variates;
        transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(control_variates), control_variate);

        double size = control_variates.size();

        double beta = 0.8;
        for(int i = 0 ; i< simulation_vector.size(); i++) {
            payoffs.push_back(m_Option.payoff(simulation_vector[i], None) - beta*discount(control_variates[i]));
        }
        double sum = 0.0, sum_of_squares = 0.0;
        for (int i =0; i<payoffs.size(); i++) {
            double discounted = discount(payoffs[i]);
            sum += discounted;
            sum_of_squares += discounted * discounted;
        }
        double estimate = size > 0 ? sum / size : NAN;
        double variance = ((sum_of_squares / size) - estimate * estimate) / size;
        return pair<double, double>(estimate+control_variate_mean, variance);
    }
}

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcDelta(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcDelta(path);
    }
    return pair<double,SensitivityMethod>(size>0? (sum/size)/sensitivityModel.deltaDivisor(): NAN, sensitivityModel.getSensitivityMethod());
};

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcGamma(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcGamma(path);
    }
    return pair<double,SensitivityMethod>(size>0? (sum/size)/sensitivityModel.gammaDivisor(): NAN, sensitivityModel.getSensitivityMethod());
};

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcVega(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcVega(path);
    }
    return pair<double,SensitivityMethod>(size>0? (sum/size)/sensitivityModel.vegaDivisor(): NAN, sensitivityModel.getSensitivityMethod());
};

#endif //SIMULATIONMETHODS_MCMODEL_H
