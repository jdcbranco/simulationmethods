
#ifndef SIMULATIONMETHODS_ASIAN_H
#define SIMULATIONMETHODS_ASIAN_H

#include "Option.h"

using namespace std;

class AsianOption: public Option {
protected:
    double m_Strike;
public:
    AsianOption(double strike, double tte, function<const double (const Path&, const Bump&)> payoff): m_Strike(strike), Option(payoff,tte) {
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
    AsianCall(double strike, double tte):
            AsianOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(path.geometric_average(bump) - m_Strike, 0.0); })
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
        return NAN;
        //The following commented code is for the arithmetic average
//
//        double P = payoff(path,None);
//        if(P>0) {
//            unsigned int m = path.size();
//            double r = params.getR();
//            double sigma = params.getSigma();
//            double sigma2 = sigma*sigma;
//            double sum = 0.0;
//            for(unsigned int i = 0; i<m; i++) {
//                double S = path.get(i);
//                double t = ((i+1)*params.getT()/m);
//                sum += S*(log(S/params.getS0()) - (r+0.5*sigma2)*t);
//            }
//            return exp(-params.getR()*params.getT()) * sum / (m*sigma);
//        } else {
//            return 0;
//        }
    };
};

/*
 * Geometric Average (G) Fixed Strike (FS) Asian call.
 */
class AsianPut: public AsianOption {
protected:
public:
    AsianPut(double strike, double tte):
            AsianOption(strike, tte, [&](const Path &path, const Bump &bump) -> double { return max(m_Strike - path.geometric_average(bump), 0.0); })
    {
    }
};

#endif //SIMULATIONMETHODS_ASIAN_H
