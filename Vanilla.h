
#ifndef SIMULATIONMETHODS_VANILLA_H
#define SIMULATIONMETHODS_VANILLA_H

#ifdef _MSC_BUILD 
#define _USE_MATH_DEFINES
#include "math.h"
inline double max(double a, double b) {
	return a > b ? a : b;
}
#endif _MSC_BUILD 

#include "Option.h"

using namespace std;

// Payoff function is passed as a lambda function.
class VanillaOption: public Option {
protected:
    double m_Strike;
public:
    VanillaOption(double strike, double tte, function<const double (const Path&, const Bump&)> payoff): m_Strike(strike), Option(payoff,tte) {
    }
    double getStrike() const {
        return m_Strike;
    }
};

class VanillaCall: public VanillaOption {
protected:
public:
    VanillaCall(double strike, double tte):
            VanillaOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(path.back(bump) - m_Strike, 0.0); }) {}
};

class VanillaPut: public VanillaOption {
protected:
public:
    VanillaPut(double strike, double tte):
            VanillaOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(m_Strike - path.back(bump), 0.0); })
    {
    }
};
#endif //SIMULATIONMETHODS_VANILLA_H
