
#ifndef SIMULATIONMETHODS_MODEL_H
#define SIMULATIONMETHODS_MODEL_H

#include <cmath>
#include "Option.h"
#include "ModelParams.h"
#include "ModelResult.h"

using namespace std;

inline double normalCDF(double value) {
    return 0.5 * erfc(-value / sqrt(2));
}

inline double normalPDF(double value) {
    return (1 / sqrt(2 * M_PI)) * exp(-0.5 * pow(value, 2));
}

/**
 * Model assumptions for now are:
 * Constant interest rate
 * Constant volatility
 * One underlying + One risk free interest rate product
 */
template<class OPTION = Option>
class Model: public ModelParams {
protected:
    OPTION &m_Option;
public:
    Model(OPTION &option, double S0, double sigma, double r): ModelParams(S0,sigma,r), m_Option(option) {}
    OPTION& getOption() const { return m_Option; }
    double getT() const override {
        return m_Option.getT();
    }
    double discount(double price) const {
        return exp(-m_r*m_Option.getT())*price;
    }

    /**
     * Calculates price and variance of the price estimate
     * @return pair where first is the price, second is the variance
     */
    virtual pair<double,double> calcPrice() const = 0;
    virtual pair<double,double> calcDelta() const = 0;
    virtual pair<double,double> calcGamma() const = 0;
    virtual pair<double,double> calcVega() const = 0;
};


#endif //SIMULATIONMETHODS_MODEL_H
