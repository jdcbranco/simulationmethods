
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
        double dt_sqrt = sqrt(dt);
        double epsilon = antithetic ? -1.0 : 1.0;
        double St = S0;
        double St_up = S0 *(1 + h); //This way, h becomes scale independent.
        double St_down = S0 * (1 - h);
        double St_sigma_up = S0;
        double St_sigma_down = S0;
        for(int i = 1; i <= random_numbers.size(); i++) {
            auto rn = epsilon * random_numbers[i-1];
            double t = i==random_numbers.size() ? T : dt*i;
            switch(model.getSolver()) {
                case Explicit: {
                    double factor = exp((r - pow(sigma, 2.0) / 2) * t + sigma * pow(t, 0.5) * rn);
                    double factor_sigma_up = exp((r - pow(sigma * (1 + h), 2.0) / 2) * t + (sigma * (1 + h)) * pow(t, 0.5) * rn);
                    double factor_sigma_down = exp((r - pow(sigma * (1 - h), 2.0) / 2) * t + (sigma * (1 - h)) * pow(t, 0.5) * rn);
                    m_Prices.push_back(S0 * factor);
                    m_Prices_bump_S_up.push_back((S0 * (1 + h)) * factor);
                    m_Prices_bump_S_down.push_back((S0 * (1 - h)) * factor);
                    m_Prices_bump_sigma_up.push_back(S0 * factor_sigma_up);
                    m_Prices_bump_sigma_down.push_back(S0 * factor_sigma_down);
                    continue;
                }
                case Euler:
                case Milstein: {
                    double St_previous = St;
                    double St_up_previous = St_up;
                    double St_down_previous = St_down;
                    double St_sigma_up_previous = St_sigma_up;
                    double St_sigma_down_previous = St_sigma_down;
                    double factor_euler = (1.0 + r * dt + sigma * dt_sqrt * rn);
                    double factor_euler_sigma_up = (1.0 + r * dt + (sigma * (1 + h)) * dt_sqrt * rn);
                    double factor_euler_sigma_down = (1.0 + r * dt + (sigma * (1 - h)) * dt_sqrt * rn);
                    St = St_previous * factor_euler;
                    St_up = St_up_previous * factor_euler;
                    St_down = St_down_previous * factor_euler;
                    St_sigma_up = St_sigma_up_previous * factor_euler_sigma_up;
                    St_sigma_down = St_sigma_down_previous * factor_euler_sigma_down;
                    m_Prices.push_back(St);
                    m_Prices_bump_S_up.push_back(St_up);
                    m_Prices_bump_S_down.push_back(St_down);
                    m_Prices_bump_sigma_up.push_back(St_sigma_up);
                    m_Prices_bump_sigma_down.push_back(St_sigma_down);
                    continue;
                }
            }
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
