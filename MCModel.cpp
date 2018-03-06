//
// Created by jdcbr on 3/4/2018.
//

#include "MCModel.h"

#include <algorithm>
#include <numeric>

using namespace std;

double MCModel::calcPrice() const {
    vector<double> payoffs;
    transform(simulation_vector.begin(), simulation_vector.end(),payoffs.begin(), m_Option.getPayoffFunction());
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