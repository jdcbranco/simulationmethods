
#ifndef SIMULATIONMETHODS_VANILLA_H
#define SIMULATIONMETHODS_VANILLA_H

#include "Option.h"

using namespace std;

// Payoff function is passed as a lambda function.
class VanillaOption: public Option {
protected:
    double m_Strike;
public:
    VanillaOption(double strike, double tte = 1.0, function<const double (const Path&)> payoff): m_Strike(strike), Option(payoff,tte) {
    }
    double getStrike() const {
        return m_Strike;
    }
};

class VanillaCall: public VanillaOption {
protected:
public:
    VanillaCall(double strike, double tte):
            VanillaOption(strike, tte, [](const Path &path) -> double { return max(path.back() - m_Strike, 0.0); })
    {
    }
};

class VanillaPut: public VanillaOption {
protected:
public:
    VanillaPut(double strike, double tte):
            VanillaOption(strike, tte, [](const Path &path) -> double { return max(m_Strike - path.back(), 0.0); })
    {
    }
};
#endif //SIMULATIONMETHODS_VANILLA_H
