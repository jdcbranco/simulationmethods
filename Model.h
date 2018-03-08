
#ifndef SIMULATIONMETHODS_MODEL_H
#define SIMULATIONMETHODS_MODEL_H

#include <cmath>
#include "Option.h"
#include "ModelParams.h"

using namespace std;


/**
 * Model assumptions for now are:
 * Constant interest rate
 * Constant volatility
 * One underlying + One risk free interest rate product
 */
class Model: public ModelParams {
protected:
    Option &m_Option;
public:
    Model(Option &option, double S0, double sigma, double r): ModelParams(S0,sigma,r), m_Option(option) {}
    Option& getOption() const { return m_Option; }
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
    virtual double calcDelta() const = 0;
    virtual double calcGamma() const = 0;
    virtual double calcVega() const = 0;
};


#endif //SIMULATIONMETHODS_MODEL_H
