//
// Created by jdcbr on 3/4/2018.
//

#ifndef SIMULATIONMETHODS_MODEL_H
#define SIMULATIONMETHODS_MODEL_H

#include <cmath>
#include "Option.h"

using namespace std;

//Model assumptions for now are:
// * Constant interest rate
// * Constant volatility
// * One underlying + risk free interest rate product
class Model {
protected:
    Option m_Option;
    double m_S0; //initial price
    double m_Sigma; //constant volatility throught the period.
    double m_r; //risk free interest rate, constant throughout the period.
public:
    Model(Option option, double S0, double sigma, double r): m_Option(option), m_S0(S0), m_Sigma(sigma), m_r(r) {}
    double discount(double price) const {
        return exp(-m_r*m_Option.getT());
    }
    Option getOption() const { return m_Option; }
    double getS0() const { return m_S0; }
    double getSigma() const { return m_Sigma; }
    double getR() const { return m_r; }
};


#endif //SIMULATIONMETHODS_MODEL_H
