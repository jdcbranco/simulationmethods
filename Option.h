
#ifndef SIMULATIONMETHODS_OPTION_H
#define SIMULATIONMETHODS_OPTION_H

#include <vector>
#include <functional>

#include "Path.h"

using namespace std;

/*
 * Base class to be used for single options.
 * The Strike and the Payoff are option-specific parameters, so they are members.
 * Interest rate and other parameters are model-specific, not option specific, thus those kind of
 * parameters don't belong here.
 * */
class Option {
protected:
    function<const double (const Path&, const Bump&)> m_Payoff;
    double m_T; //Time to Expiry
public:
    Option(function<const double (const Path&, const Bump&)> payoff, double tte): m_T(tte), m_Payoff(payoff) {}
    double payoff(const Path &path, const Bump &bump) const {
        return m_Payoff(path,bump);
    }
    double getT() const { return m_T; }
    virtual double pathwise_delta(const Path &path, const ModelParams &params) const = 0;
    virtual double pathwise_gamma(const Path &path, const ModelParams &params) const = 0;
    virtual double pathwise_vega(const Path &path, const ModelParams &params) const = 0;
};


#endif //SIMULATIONMETHODS_OPTION_H
