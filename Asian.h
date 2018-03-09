
#ifndef SIMULATIONMETHODS_ASIAN_H
#define SIMULATIONMETHODS_ASIAN_H

#include "Option.h"

using namespace std;

class AsianOption: public Option {
protected:
    double m_Strike;
public:
    AsianOption(double strike, double tte, int dim, function<const double (const Path&, const Bump&)> payoff): m_Strike(strike), Option(payoff,tte,dim) {
    }
    double getStrike() const {
        return m_Strike;
    }
};

/*
 * Geometric Average (G) Fixed Strike Asian call.
 */
class AsianCall: public AsianOption {
protected:
public:
    AsianCall(double strike, double tte, int dim):
            AsianOption(strike, tte, dim, [&](const Path &path, const Bump &bump) -> double {
                if(path.getPathType()==GeometricAverage) {
                    return max(path.back(bump) - m_Strike, 0.0);
                } else {
                    return max(path.geometric_average(bump) - m_Strike, 0.0);
                }
            })
    {
    }
};

/*
 * Geometric Average (G) Fixed Strike (FS) Asian call.
 */
class AsianPut: public AsianOption {
protected:
public:
    AsianPut(double strike, double tte, int dim):
            AsianOption(strike, tte, dim, [&](const Path &path, const Bump &bump) -> double {
                if(path.getPathType()==GeometricAverage) {
                    return max(m_Strike - path.back(bump), 0.0);
                } else {
                    return max(m_Strike - path.geometric_average(bump), 0.0);
                }})
    {
    }
};

/**
 * Arithmetic Average Fixed Strike Asian Call
 */
class ArithmeticAsianCall: public AsianOption {
protected:
public:
    ArithmeticAsianCall(double strike, double tte, int dim):
            AsianOption(strike, tte, dim, [&](const Path &path, const Bump &bump) -> double {
                if(path.getPathType()==GeometricAverage) {
                    return max(path.back(bump) - m_Strike, 0.0);
                } else {
                    return max(path.arithmetic_average(bump) - m_Strike, 0.0);
                }
            })
    {
    }
};

/**
 * Arithmetic Average Fixed Strike Asian Put
 */
class ArithmeticAsianPut: public AsianOption {
protected:
public:
    ArithmeticAsianPut(double strike, double tte, int dim):
            AsianOption(strike, tte, dim, [&](const Path &path, const Bump &bump) -> double {
                if(path.getPathType()==GeometricAverage) {
                    return max(m_Strike - path.back(bump), 0.0);
                } else {
                    return max(m_Strike - path.arithmetic_average(bump), 0.0);
                }})
    {
    }
};

#endif //SIMULATIONMETHODS_ASIAN_H
