
#ifndef SIMULATIONMETHODS_PATH_H
#define SIMULATIONMETHODS_PATH_H

#include <vector>
#include <cmath>
#include "Model.h"

using namespace std;

class Path {
protected:
    vector<double> m_Prices; //non discounted price path
public:
    Path(Model model, vector<double> &random_numbers, bool antithetic = false) {
        double S0 = model.getS0();
        double T = model.getOption().getT();
        double r = model.getR();
        double sigma = model.getSigma();
        double dt = T / random_numbers.size();
        double epsilon = antithetic ? -1.0 : 1.0;
        for(int i = 1; i <= random_numbers.size(); i++) {
            auto rn = epsilon * random_numbers[i-1];
            double t = i==random_numbers.size() ? T : dt*i;
            double S_t = S0*exp((r - pow(sigma, 2.0) / 2) * t + sigma * pow(t, 0.5) * rn);
            m_Prices.push_back(S_t);
        }
    }
};

#endif //SIMULATIONMETHODS_PATH_H
