
#include "MCModel.h"

#include <algorithm>
#include <numeric>

using namespace std;

double MCModel::calcPrice() const {
    vector<double> payoffs;
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs), [&](const Path &path) { return m_Option.payoff(path,None); });
    double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
    double size = payoffs.size();
    return size>0 ? discount(sum / size) : NAN;
}

double MCModel::calcDelta() const {
    vector<double> payoffs_bump_up;
    vector<double> payoffs_bump_down;
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up), [&](const Path &path) { return m_Option.payoff(path,Price_Up); });
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down), [&](const Path &path) { return m_Option.payoff(path,Price_Down); });
    double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
    double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
    double size = payoffs_bump_up.size();
    return size>0 ? discount((sum_up-sum_down) / size)/ (2*m_h*m_S0) : NAN; //TODO Divide before discounting or after?
}

double MCModel::calcGamma() const {
    vector<double> payoffs;
    vector<double> payoffs_bump_up;
    vector<double> payoffs_bump_down;
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs), [&](const Path &path) { return m_Option.payoff(path,None); });
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up), [&](const Path &path) { return m_Option.payoff(path,Price_Up); });
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down), [&](const Path &path) { return m_Option.payoff(path,Price_Down); });
    double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
    double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
    double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
    double size = payoffs_bump_up.size();
    return size>0 ? discount((sum_up+sum_down-2.0*sum) / size)/ (m_h*m_h*m_S0*m_S0) : NAN;
}

double MCModel::calcVega() const {
    vector<double> payoffs_bump_up;
    vector<double> payoffs_bump_down;
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_up), [&](const Path &path) { return m_Option.payoff(path,Sigma_Up); });
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs_bump_down), [&](const Path &path) { return m_Option.payoff(path,Sigma_Down); });
    double sum_up = accumulate(payoffs_bump_up.begin(), payoffs_bump_up.end(), 0.0);
    double sum_down = accumulate(payoffs_bump_down.begin(), payoffs_bump_down.end(), 0.0);
    double size = payoffs_bump_up.size();
    return size>0 ? discount((sum_up-sum_down) / size)/ (2*m_h*m_Sigma) : NAN;
}