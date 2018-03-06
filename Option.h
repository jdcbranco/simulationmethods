//
// Created by jdcbr on 3/4/2018.
//

#ifndef SIMULATIONMETHODS_OPTION_H
#define SIMULATIONMETHODS_OPTION_H

#include <vector>
#include <functional>

#include "Path.h"

using namespace std;

//This is necessary to to avoid circular reference, thus letting Option reference Path before Path is actually defined.
//class Path;

class Option {
protected:
    function<const double (const Path&)> m_Payoff;
    double m_T; //Time to Expiry
public:
    Option(function<const double (const Path&)> payoff, double tte): m_T(tte), m_Payoff(payoff) {}
    double payoff(const Path &path) const {
        return m_Payoff(path);
    }
    double getT() const { return m_T; }
    function<const double(const Path&)> getPayoffFunction() const { return m_Payoff; }
};


#endif //SIMULATIONMETHODS_OPTION_H
