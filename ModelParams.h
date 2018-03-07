
#ifndef SIMULATIONMETHODS_MODELPARAMS_H
#define SIMULATIONMETHODS_MODELPARAMS_H

using namespace std;

enum SDESolver { Explicit, Euler, Milstein };

class ModelParams {
protected:
    double m_S0; //initial price
    double m_Sigma; //constant volatility throught the period.
    double m_r; //risk free interest rate, constant throughout the period.
    double m_h; //bump in Prices and Sigma to calculate the greeks
    SDESolver m_Solver;
public:
    ModelParams(double S0, double sigma, double r, double h = 0.01, SDESolver sdeSolver = Explicit): m_S0(S0), m_Sigma(sigma), m_r(r), m_h(h), m_Solver(sdeSolver) {}
    virtual double getT() const = 0;
    double getS0() const { return m_S0; }
    double getSigma() const { return m_Sigma; }
    double getR() const { return m_r; }
    double getH() const { return m_h; }
    SDESolver getSolver() const { return m_Solver; }
};


#endif //SIMULATIONMETHODS_MODELPARAMS_H
