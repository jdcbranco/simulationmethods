
#ifndef SIMULATIONMETHODS_PATH_H
#define SIMULATIONMETHODS_PATH_H

#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include "ModelParams.h"

using namespace std;
static bool print_sample_path = false;
enum Bump { None, Price_Up, Price_Down, Sigma_Up, Sigma_Down };
enum PathType { Price, GeometricAverage };

/**
 * Stores the price paths and any pertubation to them required for finite difference methods.
 * By perturbing the path using the same random numbers, we save computation time.
 * We can argue that this means more memory needs, but we can perform optimizations that
 * discard and free the path once all calculations to them have been done, so we could in theory
 * process the Paths in batches, either serially or in parallel.
 * Note that this code is highly couple to the Model code - assumes constant interest rate and volatility.
 * If we want to make it totally flexible, we must inject here the drift and diffusion equation formulas,
 * and explicit solution formulas, if they exist.
 */
class Path {
private:
    //TODO It doesn't feel right to let this function be here. May move it to some utils class
    double geometric_average(vector<double> input) const {
        vector<double> logPrices;
        double acc = 0.0;
        for(auto item: input) {
            acc += log(item);
        }
        return exp(acc / input.size());
    }
    double arithmetic_average(vector<double> input) const {
        return input.size()>0 ? accumulate(input.begin(), input.end(), 0.0) / input.size() : NAN;
    }
protected:
    PathType m_PathType;
    double m_S0, m_Epsilon;
    vector<double> m_Prices; //non discounted price path
    vector<double> m_Prices_bump_S_up;
    vector<double> m_Prices_bump_S_down;
    vector<double> m_Prices_bump_sigma_up;
    vector<double> m_Prices_bump_sigma_down;
    vector<double> m_RandomNumbers;
public:
    Path(ModelParams &model, vector<double> &&random_numbers, bool antithetic = false):
            m_RandomNumbers(random_numbers),
            m_S0(model.getS0()),
            m_PathType(model.getSolver()==ExplicitGeometricAverage?GeometricAverage:Price) {
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
        m_Epsilon = antithetic ? -1.0 : 1.0;
        //For use in the Euler method
        double St = S0;
        double St_up = S0*(1+h); //This way, h becomes scale independent.
        double St_down = S0*(1-h);
        double St_sigma_up = S0;
        double St_sigma_down = S0;
        //For use in the Explicit method
        double factor_explicit = 0.0, factor_explicit_sigma_up = 0.0, factor_explicit_sigma_down = 0.0;
        for(int i = 1; i <= random_numbers.size(); i++) {
            auto rn = m_Epsilon * random_numbers[i-1];
            switch(model.getSolver()) {
                case Explicit: {
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
                case ExplicitGeometricAverage: {
                    double n = model.getDim();
                    double a = (n+1)/(2*n);
                    double b = sqrt((n+1)*(2*n+1)/(6*n*n));
                    factor_explicit += (r - sigma2 / 2) * a * dt + sigma * b * dt_sqrt * rn;
                    factor_explicit_sigma_up += (r - sigma2up / 2) * a * dt +  sigma_up * b * dt_sqrt * rn;
                    factor_explicit_sigma_down += (r - sigma2down / 2) * a * dt +  sigma_down * b * dt_sqrt * rn;
                    m_Prices.push_back(S0 * exp(factor_explicit));
                    m_Prices_bump_S_up.push_back((S0 * (1 + h)) * exp(factor_explicit));
                    m_Prices_bump_S_down.push_back((S0 * (1 - h)) * exp(factor_explicit));
                    m_Prices_bump_sigma_up.push_back(S0 * exp(factor_explicit_sigma_up));
                    m_Prices_bump_sigma_down.push_back(S0 * exp(factor_explicit_sigma_down));
                    continue;
                }
                case Euler: {
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
    }

    PathType getPathType() const {
        return m_PathType;
    }
    unsigned int size() const {
        return m_Prices.size();
    }
    double get(unsigned int i) const {
        return m_Prices[i];
    }
    double front_random_number() const {
        return m_RandomNumbers.size()>0 ? m_Epsilon * m_RandomNumbers.front() : NAN;
    }
    double random_number(unsigned int i) const {
        return m_RandomNumbers.size()>i ? m_Epsilon * m_RandomNumbers[i] : NAN;
    }
    double back() const {
        return back(None);
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
    double arithmetic_average(Bump bump) const {
        switch (bump) {
            case None:
                return arithmetic_average(m_Prices);
            case Price_Up:
                return arithmetic_average(m_Prices_bump_S_up);
            case Price_Down:
                return arithmetic_average(m_Prices_bump_S_down);
            case Sigma_Up:
                return arithmetic_average(m_Prices_bump_sigma_up);
            case Sigma_Down:
                return arithmetic_average(m_Prices_bump_sigma_down);
        }
    }
    friend ostream& operator<<(ostream& os, const Path &path) {
        cout << "-------------" << endl;
        cout << "Sample path: " << endl;
        cout << path.m_S0 << ", ";
        for (int i = 0; i < path.m_Prices.size(); ++i) {
            cout << path.m_Prices[i] << ", " ;
        }
        cout << endl;
        return os;
    };
};

#endif //SIMULATIONMETHODS_PATH_H
