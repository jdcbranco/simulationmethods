
#ifndef SIMULATIONMETHODS_PATH_H
#define SIMULATIONMETHODS_PATH_H

#include <vector>
#include <cmath>
#include "ModelParams.h"

using namespace std;

enum Bump { None, Price_Up, Price_Down, Sigma_Up, Sigma_Down };

class Path {
protected:
    vector<double> m_Prices; //non discounted price path
    vector<double> m_Prices_bump_S_up;
    vector<double> m_Prices_bump_S_down;
    vector<double> m_Prices_bump_sigma_up;
    vector<double> m_Prices_bump_sigma_down;
    vector<double> &m_RandomNumbers;
public:
    Path(ModelParams &model, vector<double> &&random_numbers, bool antithetic = false): m_RandomNumbers(random_numbers) {
        double S0 = model.getS0();
        double T = model.getT();
        double r = model.getR();
        double h = model.getH();//Bump in the parameters in order to calculate the sensitivities
        double sigma = model.getSigma();
        double dt = T / random_numbers.size();
        double epsilon = antithetic ? -1.0 : 1.0;
        for(int i = 1; i <= random_numbers.size(); i++) {
            auto rn = epsilon * random_numbers[i-1];
            double t = i==random_numbers.size() ? T : dt*i;
            double factor = exp((r - pow(sigma, 2.0) / 2) * t + sigma * pow(t, 0.5) * rn);
            double factor_sigma_up = exp((r - pow(sigma+h, 2.0) / 2) * t + (sigma+h) * pow(t, 0.5) * rn);
            double factor_sigma_down = exp((r - pow(sigma-h, 2.0) / 2) * t + (sigma-h) * pow(t, 0.5) * rn);

            double S_t = S0*factor;
            m_Prices.push_back(S0*factor);
            m_Prices_bump_S_up.push_back((S0+h)*factor);
            m_Prices_bump_S_down.push_back((S0-h)*factor);
            m_Prices_bump_sigma_up.push_back(S0*factor_sigma_up);
            m_Prices_bump_sigma_down.push_back(S0*factor_sigma_down);
        }
    }
    double back(Bump bump) const {
        switch (bump) {
            case None:
                return m_Prices.back();
            case Price_Up:
                return m_Prices_bump_S_up.back();
            case Price_Down:
                return m_Prices_bump_S_down.back();
            case Sigma_Up:
                return m_Prices_bump_sigma_up.back();
            case Sigma_Down:
                return m_Prices_bump_sigma_down.back();
        }
    }
};

#endif //SIMULATIONMETHODS_PATH_H
