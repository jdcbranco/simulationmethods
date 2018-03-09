
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

    ModelResult simulate(Simulator simulator, int simulations, int path_size = 1, SensitivityMethod sensitivityMethod = SensitivityMethod::FiniteDifference) {
        clock_t start = clock();
        this->simulation_vector.clear();
        this->simulation_vector = simulator.simulate(*this, simulations, path_size);
        this->m_SensitivityMethod = sensitivityMethod;
        auto price_and_variance = this->calcPrice();
        auto delta = this->calcDelta();
        auto gamma = this->calcGamma();
        auto vega  = this->calcVega();
        ModelResult result;
        result.setModelType(ModelType::MonteCarlo);
        result.setDeltaMethod(delta.second);
        result.setGammaMethod(gamma.second);
        result.setVegaMethod(vega.second);
        result.setPrice(price_and_variance.first);
        result.setPriceVariance(price_and_variance.second);
        result.setDelta(delta.first);
        result.setGamma(gamma.first);
        result.setVega(vega.first);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
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
    pair<double,SensitivityMethod> calcDelta() const override ;
    pair<double,SensitivityMethod> calcGamma() const override;
    pair<double,SensitivityMethod> calcVega() const override;
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
pair<double,SensitivityMethod> MCModel<OPTION>::calcDelta() const {
    switch (m_SensitivityMethod) {
        case SensitivityMethod::FiniteDifference: {
            vector<double> payoffs_bump_up;
            vector<double> payoffs_bump_down;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up),
                      [&](const Path &path) { return m_Option.payoff(path, Price_Up); });
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down),
                      [&](const Path &path) { return m_Option.payoff(path, Price_Down); });
            double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
            double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
            double size = payoffs_bump_up.size();
            return pair<double,SensitivityMethod>(size > 0 ? discount((sum_up - sum_down) / size) / (2 * m_h * m_S0) : NAN, SensitivityMethod::FiniteDifference);
        }
        case SensitivityMethod::PathwiseDifferentiation: {
            vector<double> pathwise_deltas;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(pathwise_deltas),
                      [&](const Path &path) {
                          return m_Option.pathwise_delta(path,*this);
                      });
            double sum = accumulate(pathwise_deltas.begin(), pathwise_deltas.end(), 0.0);
            double size = pathwise_deltas.size();
            return pair<double,SensitivityMethod>(size > 0 ? sum / size : NAN, SensitivityMethod::PathwiseDifferentiation);
        }
        case SensitivityMethod::LikelihoodRatio: {
            vector<double> payoffs;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs),
                      [&](const Path &path) {
                          double S_T = path.back(None);
                          double discounted_payoff = discount(m_Option.payoff(path, None));
                          double dt_sqrt = sqrt(m_Option.getT()/path.size());
                          double Z = 0;
                          for(int i=0; i<path.size(); i++) {
                              Z += path.random_number(i);
                          }
                          Z /= sqrt(path.size());
                          return discounted_payoff * Z / (m_S0 * m_Sigma * sqrt(m_Option.getT()));
                      });
            double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
            double size = payoffs.size();
            return pair<double,SensitivityMethod>(size > 0? sum / size : NAN, SensitivityMethod::LikelihoodRatio);
        }
    }
}

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcGamma() const {
    switch(m_SensitivityMethod) {
        case SensitivityMethod::FiniteDifference: {
            vector<double> payoffs;
            vector<double> payoffs_bump_up;
            vector<double> payoffs_bump_down;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs),
                      [&](const Path &path) { return m_Option.payoff(path, None); });
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up),
                      [&](const Path &path) { return m_Option.payoff(path, Price_Up); });
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down),
                      [&](const Path &path) { return m_Option.payoff(path, Price_Down); });
            double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
            double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
            double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
            double size = payoffs_bump_up.size();
            return pair<double,SensitivityMethod>(size > 0 ? discount((sum_up + sum_down - 2.0 * sum) / size) / (m_h * m_h * m_S0 * m_S0) : NAN, SensitivityMethod::FiniteDifference);
        }
        case SensitivityMethod::PathwiseDifferentiation: {
            vector<double> pathwise_gamma;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(pathwise_gamma),
                      [&](const Path &path) {
                          return m_Option.pathwise_gamma(path,*this);
                      });
            double sum = accumulate(pathwise_gamma.begin(), pathwise_gamma.end(), 0.0);
            double size = pathwise_gamma.size();
            return pair<double,SensitivityMethod>(size > 0 ? sum / size : NAN, SensitivityMethod::PathwiseDifferentiation);
        }
        case SensitivityMethod::LikelihoodRatio: {
            vector<double> payoffs;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs),
                      [&](const Path &path) {
                          double S_T = path.back(None);
                          double discounted_payoff = discount(m_Option.payoff(path, None));
                          double Z = path.front_random_number();
                          return discounted_payoff * ((Z*Z-1)/(m_S0*m_S0*m_Sigma*m_Sigma*m_Option.getT()) - Z / (m_S0 * m_S0* m_Sigma * sqrt(m_Option.getT())));
                      });
            double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
            double size = payoffs.size();
            return pair<double,SensitivityMethod>(size > 0? sum / size : NAN, SensitivityMethod::LikelihoodRatio);
        }

    }
}

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcVega() const {
    switch (m_SensitivityMethod) {
        case SensitivityMethod::FiniteDifference: {
            vector<double> payoffs_bump_up;
            vector<double> payoffs_bump_down;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up),
                      [&](const Path &path) { return m_Option.payoff(path, Sigma_Up); });
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down),
                      [&](const Path &path) { return m_Option.payoff(path, Sigma_Down); });
            double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
            double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
            double size = payoffs_bump_up.size();
            return pair<double,SensitivityMethod>(size > 0 ? discount((sum_up - sum_down) / size) / (2 * m_h * m_Sigma) : NAN, SensitivityMethod::FiniteDifference);
        }
        case SensitivityMethod::PathwiseDifferentiation: {
            vector<double> pathwise_vega;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(pathwise_vega),
                      [&](const Path &path) {
                          return m_Option.pathwise_vega(path,*this);
                      });
            double sum = accumulate(pathwise_vega.begin(), pathwise_vega.end(), 0.0);
            double size = pathwise_vega.size();
            return pair<double,SensitivityMethod>(size > 0 ? sum / size : NAN, SensitivityMethod::PathwiseDifferentiation);
        }
        case SensitivityMethod::LikelihoodRatio: {
            vector<double> payoffs;
            transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs),
                      [&](const Path &path) {
                          double S_T = path.back(None);
                          double discounted_payoff = discount(m_Option.payoff(path, None));
                          double score = 0.0;
                          double dt_sqrt = sqrt(m_Option.getT()/path.size());
                          for(int i=0; i<path.size(); i++) {
                              double Z = path.front_random_number();
                              score += (Z*Z-1.0)/m_Sigma - Z * dt_sqrt;
                          }

                          return discounted_payoff * score;
                      });

            double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
            double size = payoffs.size();
            return pair<double,SensitivityMethod>(size > 0? sum / size : NAN, SensitivityMethod::LikelihoodRatio);
        }
    }
}

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcDelta(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcDelta(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcGamma(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcGamma(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};

template<class OPTION>
pair<double,SensitivityMethod> MCModel<OPTION>::calcVega(SensitivityModel<OPTION> &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcVega(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};

#endif //SIMULATIONMETHODS_MCMODEL_H
