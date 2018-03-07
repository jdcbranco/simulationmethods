
#ifndef SIMULATIONMETHODS_ASIAN_H
#define SIMULATIONMETHODS_ASIAN_H

#include "Option.h"

using namespace std;

class FixedStrikeAsianOption: public Option {
protected:
    double m_Strike;
public:
    FixedStrikeAsianOption(double strike, double tte, function<const double (const Path&, const Bump&)> payoff): m_Strike(strike), Option(payoff,tte) {
    }
    double getStrike() const {
        return m_Strike;
    }
};

/*
 * Geometric Average (G) Fixed Strike Asian call.
 */
class GFixedStrikeAsianCall: public FixedStrikeAsianOption {
protected:
public:
    GFixedStrikeAsianCall(double strike, double tte):
            FixedStrikeAsianOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(path.geometric_average(bump) - m_Strike, 0.0); })
    {
    }
};

/*
 * Geometric Average (G) Fixed Strike (FS) Asian call.
 */
class GFixedStrikeAsianPut: public FixedStrikeAsianOption {
protected:
public:
    GFixedStrikeAsianPut(double strike, double tte):
            FixedStrikeAsianOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(m_Strike - path.geometric_average(bump), 0.0); })
    {
    }
};

#endif //SIMULATIONMETHODS_ASIAN_H
