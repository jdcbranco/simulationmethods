
#ifndef SIMULATIONMETHODS_PATH_H
#define SIMULATIONMETHODS_PATH_H

#include <vector>
#include <cmath>
#include <iostream>
#include "ModelParams.h"

using namespace std;
static bool print_sample_path = false;
enum Bump { None, Price_Up, Price_Down, Sigma_Up, Sigma_Down };

class Path {
private:
    double geometric_average(vector<double> input) const {
        vector<double> logPrices;
        double acc = 0.0;
        for(auto item: input) {
            acc += log(item);
        }
        return exp(acc / input.size());
    }
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
        double sigma_up = sigma*(1+h);
        double sigma_down = sigma*(1-h);
        double sigma2 = sigma*sigma;
        double sigma2up = sigma_up*sigma_up;
        double sigma2down = sigma_down*sigma_down;
        double dt = T / random_numbers.size();
        double dt_sqrt = sqrt(dt);
        double epsilon = antithetic ? -1.0 : 1.0;
        //For use in the Euler method
        double St = S0;
        double St_up = S0*(1+h); //This way, h becomes scale independent.
        double St_down = S0*(1-h);
        double St_sigma_up = S0;
        double St_sigma_down = S0;
        //For use in the Explicit method
        double factor_explicit = 0.0, factor_explicit_sigma_up = 0.0, factor_explicit_sigma_down = 0.0;
        for(int i = 1; i <= random_numbers.size(); i++) {
            auto rn = epsilon * random_numbers[i-1];
            double t = i==random_numbers.size() ? T : dt*i;
            double t_sqrt = sqrt(t);
            switch(model.getSolver()) {
                case Explicit:
                {
                    factor_explicit += (r - sigma2 / 2) * dt + sigma * dt_sqrt * rn;
                    factor_explicit_sigma_up += (r - sigma2up / 2) * dt +  sigma_up * dt_sqrt * rn;
                    factor_explicit_sigma_down += (r - sigma2down / 2) * dt +  sigma_down * dt_sqrt * rn;
                    m_Prices.push_back(S0 * exp(factor_explicit));
                    m_Prices_bump_S_up.push_back((S0 * (1 + h)) * exp(factor_explicit));
                    m_Prices_bump_S_down.push_back((S0 * (1 - h)) * exp(factor_explicit));
                    m_Prices_bump_sigma_up.push_back(S0 * exp(factor_explicit_sigma_up));
                    m_Prices_bump_sigma_down.push_back(S0 * exp(factor_explicit_sigma_down));
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
                    double factor_euler_sigma_up = (1.0 + r * dt + sigma_up * dt_sqrt * rn);
                    double factor_euler_sigma_down = (1.0 + r * dt + sigma_down * dt_sqrt * rn);
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
        if(print_sample_path) {
            cout << "-------------" << endl;
            cout << "Sample path: " << endl;
            cout << S0 << ", ";
            for (int i = 0; i < m_Prices.size(); ++i) {
                cout << m_Prices[i] << ", " ;
            }
            cout << endl;
            print_sample_path=false;
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
    double geometric_average(Bump bump) const {
        switch (bump) {
            case None:
                return geometric_average(m_Prices);
            case Price_Up:
                return geometric_average(m_Prices_bump_S_up);
            case Price_Down:
                return geometric_average(m_Prices_bump_S_down);
            case Sigma_Up:
                return geometric_average(m_Prices_bump_sigma_up);
            case Sigma_Down:
                return geometric_average(m_Prices_bump_sigma_down);
        }
    }
};

#endif //SIMULATIONMETHODS_PATH_H
