
#ifndef SIMULATIONMETHODS_BSASIAN_H
#define SIMULATIONMETHODS_BSASIAN_H

#include <ctime>
#include "Model.h"
#include "Asian.h"
#include "ModelResult.h"

class BSAsianModel: public Model {
protected:
    double m_K;
public:
    BSAsianModel(FixedStrikeAsianOption &option, double S0, double sigma, double r): Model(option,S0,sigma,r), m_K(option.getStrike()) {}
    ModelResult calculate() {
        clock_t start = clock();
        auto price_and_variance = this->calcPrice();
        auto delta = this->calcDelta();
        auto gamma = this->calcGamma();
        auto vega  = this->calcVega();
        ModelResult result;
        result.setModelType(ModelType::Analytical);
        result.setDeltaMethod(delta.second);
        result.setGammaMethod(gamma.second);
        result.setVegaMethod(vega.second);
        result.setPrice(price_and_variance.first);
        result.setPriceVariance(price_and_variance.second);
        result.setDelta(delta.first);
        result.setGamma(gamma.first);
        result.setVega(vega.first);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
};

/**
 * Black Scholes model (based on BS model assumptions) for Geometric Average Fixed Strike Asian call options.
 */
class BSAsianCallModel: public BSAsianModel {
public:
    BSAsianCallModel(GFixedStrikeAsianCall &call, double S0, double sigma, double r, double path_size): BSAsianModel(call,S0,sigma,r) {}
    pair<double,double> calcPrice() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        double d2 = d1 - m_Sigma * sqrt(m_Option.getT());
        return pair<double,double>(m_S0*normalCDF(d1) - normalCDF(d2)*m_K*exp(-m_r*T),0.0);
    };
    pair<double,SensitivityMethod> calcDelta() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,SensitivityMethod>(normalCDF(d1),SensitivityMethod::Analytical);
    };
    pair<double,SensitivityMethod> calcGamma() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,SensitivityMethod>(exp(-d1*d1 / 2) / sqrt(2 * M_PI) / (m_S0*m_Sigma*sqrt(T)), SensitivityMethod::Analytical);
    };
    pair<double,SensitivityMethod> calcVega() const override {
        double T = m_Option.getT();
        double d1 = (log(m_S0 / m_K) + (m_r + m_Sigma*m_Sigma / 2)*T) / (m_Sigma*sqrt(T));
        return pair<double,SensitivityMethod>(exp(-d1*d1 / 2) / sqrt(2 * M_PI) * m_S0 *sqrt(T), SensitivityMethod::Analytical);
    };
};

#endif //SIMULATIONMETHODS_BSASIAN_H
