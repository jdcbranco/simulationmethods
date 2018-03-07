
#ifndef SIMULATIONMETHODS_BSVANILLA_H
#define SIMULATIONMETHODS_BSVANILLA_H

#include <ctime>
#include "Model.h"
#include "Vanilla.h"
#include "ModelResult.h"

using namespace std;

double normalCDF(double value) {
    return 0.5 * erfc(-value / sqrt(2));
}

class BSModel: public Model {
protected:
    double m_K;
public:
    BSModel(VanillaOption &option, double S0, double sigma, double r): Model(option,S0,sigma,r), m_K(option.getStrike()) {}
    ModelResult calculate() {
        clock_t start = clock();
        double price = this->calcPrice();
        double delta = this->calcDelta();
        double gamma = this->calcGamma();
        double vega  = this->calcVega();
        ModelResult result;
        result.setPrice(price);
        result.setDelta(delta);
        result.setGamma(gamma);
        result.setVega(vega);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
};

/**
 * Black Scholes model for European call options.
 * Based on original version from european.cpp
 */
class BSCallModel: public BSModel {
public:
    BSCallModel(VanillaCall &call, double S0, double sigma, double r): BSModel(call,S0,sigma,r) {}
    double calcPrice() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        double d2 = d1 - m_Sigma * sqrt(m_Option.getT());
        return (m_S0*normalCDF(d1) - normalCDF(d2)*m_K*exp(-m_r*T));
    };
    double calcDelta() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return (normalCDF(d1));
    };
    double calcGamma() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return (exp(-d1*d1 / 2) / sqrt(2 * M_PI) / (m_S0*m_Sigma*sqrt(T)));
    };
    double calcVega() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return(exp(-d1*d1 / 2) / sqrt(2 * M_PI) * m_S0 *sqrt(T));
    };
};

/**
 * Black Scholes model for European put options.
 */
//TODO Implement this class
class BSPutModel: public BSModel {
public:
    BSPutModel(VanillaPut &put, double S0, double sigma, double r): BSModel(put,S0,sigma,r) {}
};


#endif //SIMULATIONMETHODS_BSVANILLA_H
