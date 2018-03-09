
#include "MCModel.h"

#include <algorithm>
#include <numeric>

using namespace std;

pair<double,double> MCModel::calcPrice() const {
    if(!control_variate) {
        vector<double> payoffs;
        transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs),
                  [&](const Path &path) { return m_Option.payoff(path, None); });
        double sum = 0.0, sum_of_squares = 0.0;
        for (auto i = payoffs.begin(); i != payoffs.end(); i++) {
            double discounted = discount(*i);
            sum += discounted;
            sum_of_squares += discounted * discounted;
        }
        double size = payoffs.size();
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

pair<double,SensitivityMethod> MCModel::calcDelta() const {
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

pair<double,SensitivityMethod> MCModel::calcGamma() const {
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

pair<double,SensitivityMethod> MCModel::calcVega() const {
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

pair<double,SensitivityMethod> MCModel::calcDelta(SensitivityModel &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcDelta(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};
pair<double,SensitivityMethod> MCModel::calcGamma(SensitivityModel &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcGamma(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};
pair<double,SensitivityMethod> MCModel::calcVega(SensitivityModel &sensitivityModel) const {
    double sum = 0.0;
    double size = simulation_vector.size();
    for (auto path : simulation_vector) {
        sum += sensitivityModel.calcVega(path);
    }
    return pair<double,SensitivityMethod>(size>0? sum/size: NAN, sensitivityModel.getSensitivityMethod());
};