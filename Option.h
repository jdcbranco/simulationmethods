
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
    int m_dim;
public:
    Option(function<const double (const Path&, const Bump&)> payoff, double tte, int dim = 1): m_T(tte), m_Payoff(payoff), m_dim(dim) {}
    double payoff(const Path &path) const {
        return payoff(path, None);
    }
    double payoff(const Path &path, const Bump &bump) const {
        return m_Payoff(path,bump);
    }
    double getT() const { return m_T; }
    int getDim() const { return m_dim; }
};


#endif //SIMULATIONMETHODS_OPTION_H
