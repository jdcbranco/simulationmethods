
#ifndef SIMULATIONMETHODS_BSVANILLA_H
#define SIMULATIONMETHODS_BSVANILLA_H

#include <ctime>
#include "Model.h"
#include "Vanilla.h"
#include "ModelResult.h"

using namespace std;

template<class OPTION = VanillaOption>
class BSModel: public Model<OPTION> {
protected:
    double m_K;
public:
    BSModel(OPTION &option, double S0, double sigma, double r): Model<OPTION>(option,S0,sigma,r), m_K(option.getStrike()) {}
    ModelResult calculate() {
        clock_t start = clock();
        auto price = this->calcPrice();
        auto delta = this->calcDelta();
        auto gamma = this->calcGamma();
        auto vega  = this->calcVega();
        ModelResult result;
        result.setModelType(ModelType::Analytical);
        result.setAntitheticVariate(false);
        result.setControlVariate(false);
        result.setGreeksMethod(SensitivityMethod::Analytical);
        result.setPrice(price.first);
        result.setPriceVariance(price.second);
        result.setDelta(delta.first);
        result.setDeltaVariance(delta.second);
        result.setGamma(gamma.first);
        result.setGammaVariance(gamma.second);
        result.setVega(vega.first);
        result.setVegaVariance(vega.second);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
};

/**
 * Black Scholes model for European call options.
 * Based on original version from european.cpp
 */
class BSCallModel: public BSModel<VanillaCall> {
public:
    BSCallModel(VanillaCall &call, double S0, double sigma, double r): BSModel(call,S0,sigma,r) {}
    pair<double,double> calcPrice() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        double d2 = d1 - m_Sigma * sqrt(m_Option.getT());
        return pair<double,double>(m_S0*normalCDF(d1) - normalCDF(d2)*m_K*exp(-m_r*T),0.0);
    };
    pair<double,double> calcDelta() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,double>(normalCDF(d1), 0.0);
    };
    pair<double,double> calcGamma() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,double>(exp(-d1*d1 / 2) / sqrt(2 * M_PI) / (m_S0*m_Sigma*sqrt(T)), 0.0);
    };
    pair<double,double> calcVega() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,double>(exp(-d1*d1 / 2) / sqrt(2 * M_PI) * m_S0 *sqrt(T), 0.0);
    };
};

/**
 * Black Scholes model for European put options.
 */
//TODO Implement this class
class BSPutModel: public BSModel<VanillaPut> {
public:
    BSPutModel(VanillaPut &put, double S0, double sigma, double r): BSModel(put,S0,sigma,r) {}
};


#endif //SIMULATIONMETHODS_BSVANILLA_H
