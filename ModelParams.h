
#ifndef SIMULATIONMETHODS_MODELPARAMS_H
#define SIMULATIONMETHODS_MODELPARAMS_H

using namespace std;

class ModelParams {
protected:
    double m_S0; //initial price
    double m_Sigma; //constant volatility throught the period.
    double m_r; //risk free interest rate, constant throughout the period.
public:
    ModelParams(double S0, double sigma, double r): m_S0(S0), m_Sigma(sigma), m_r(r) {}
    virtual double getT() const = 0;
    double getS0() const { return m_S0; }
    double getSigma() const { return m_Sigma; }
    double getR() const { return m_r; }
};


#endif //SIMULATIONMETHODS_MODELPARAMS_H
