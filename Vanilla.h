
#ifndef SIMULATIONMETHODS_VANILLA_H
#define SIMULATIONMETHODS_VANILLA_H

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
            VanillaOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(path.back(bump) - m_Strike, 0.0); })
    {
    }
    double pathwise_delta(const Path &path, const ModelParams &params) const override {
        double P = payoff(path,None);
        return exp(-params.getR()*params.getT())* (P>0? (P+m_Strike)/params.getS0(): 0.0);
    };
    double pathwise_gamma(const Path &path, const ModelParams &params) const override {
        return NAN;
    };
    double pathwise_vega(const Path &path, const ModelParams &params) const override {
        double P = payoff(path,None);
        if(P>0) {
            double r = params.getR();
            double S = P+m_Strike;
            double sigma = params.getSigma();
            double sigma2 = sigma*sigma;
            return exp(-params.getR()*params.getT())* (S/sigma)*(log(S/params.getS0()) - (r+0.5*sigma2)*params.getT());
        } else {
            return 0;
        }
    };
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
