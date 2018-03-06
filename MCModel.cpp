//
// Created by jdcbr on 3/4/2018.
//

#include "MCModel.h"

#include <algorithm>
#include <numeric>

using namespace std;

double MCModel::calcPrice() const {
    vector<double> payoffs;
    //m_Option.getPayoffFunction()
    transform(simulation_vector.begin(), simulation_vector.end(), back_inserter(payoffs), [&](const Path &path) { return m_Option.payoff(path); });
    double sum = accumulate(payoffs.begin(), payoffs.end(), 0.0);
    double size = payoffs.size();
    return size>0 ? discount(sum / payoffs.size()) : NAN;
}

double MCModel::calcDelta() const {

}

double MCModel::calcGamma() const {

}

double MCModel::calcVega() const {

}