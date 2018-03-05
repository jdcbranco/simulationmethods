//
// Created by jdcbr on 3/4/2018.
//

#ifndef SIMULATIONMETHODS_OPTION_H
#define SIMULATIONMETHODS_OPTION_H

#include <vector>
#include <functional>

using namespace std;

class Option {
protected:
    function<double (const vector<double>&)> m_Payoff;
    double m_T; //Time to Expiry
public:
    Option(function<double (const vector<double>&)> payoff, double tte): m_T(tte), m_Payoff(payoff) {}
    double payoff(vector<double> path) const {
        return m_Payoff(path);
    }
    double getT() const { return m_T; }
};


#endif //SIMULATIONMETHODS_OPTION_H
