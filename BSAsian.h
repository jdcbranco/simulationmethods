
#ifndef SIMULATIONMETHODS_BSASIAN_H
#define SIMULATIONMETHODS_BSASIAN_H

#include <ctime>
#include "Model.h"
#include "Asian.h"
#include "ModelResult.h"

template<class OPTION = AsianOption>
class BSAsianModel: public Model<OPTION> {
protected:
    double m_K;
    double m_PathSize;
public:
    BSAsianModel(OPTION &option, double S0, double sigma, double r, double path_size): Model<OPTION>(option,S0,sigma,r), m_K(option.getStrike()), m_PathSize(path_size) {}
    ModelResult calculate() {
        clock_t start = clock();
        auto price = this->calcPrice();
        auto delta = this->calcDelta();
        auto gamma = this->calcGamma();
        auto vega  = this->calcVega();
        ModelResult result;
        result.setModelType(ModelType::Analytical);
        result.setAntitheticVariate(false);
        result.setControlVariate(false);
        result.setGreeksMethod(SensitivityMethod::Analytical);
        result.setPrice(price.first);
        result.setPriceVariance(price.second);
        result.setDelta(delta.first);
        result.setDeltaVariance(delta.second);
        result.setGamma(gamma.first);
        result.setGammaVariance(gamma.second);
        result.setVega(vega.first);
        result.setVegaVariance(vega.second);
        result.setCalcTime((std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
        return result;
    }
};

/**
 * Black Scholes model (based on BS model assumptions) for Geometric Average Fixed Strike Asian call options.
 */
class BSAsianCallModel: public BSAsianModel<AsianCall> {
public:
    BSAsianCallModel(AsianCall &call, double S0, double sigma, double r, double path_size): BSAsianModel(call,S0,sigma,r,path_size) {}
    pair<double,double> calcPrice() const override {
        double T = m_Option.getT();
        double v = m_Sigma;
        double r = m_r;
        double N = m_PathSize;
        double K = m_K;
        double S0 = m_S0;

        double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;
        double discount  = exp(-r*T);
        double price = discount*(S0 * exp(mu_a + 0.5*var_a)*normalCDF(d1) - K*normalCDF(d2));

        return pair<double,double>(price,0.0);
    };
    pair<double,double> calcDelta() const override {
        double T = m_Option.getT();
        double r = m_r;
        double v = m_Sigma;
        double t = 0;
        double S0 = m_S0;
        double K = m_K;
        double N = m_PathSize;

        double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;
        double delta = discount(exp(mu_a+0.5*var_a)*(normalCDF(d1)+normalPDF(d1)/sd_a) - K*normalPDF(d2)/(sd_a*S0));

        return pair<double,double>(delta, 0.0);
    };
    pair<double,double> calcGamma() const override {
        double T = m_Option.getT();
        double r = m_r;
        double v = m_Sigma;
        double t = 0;
        double S0 = m_S0;
        double K = m_K;
        double N = m_PathSize;

        double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;
        double d_d = 1.0 / (sd_a*S0);
        double gamma =  discount(exp( mu_a + pow(sd_a,2 ) ) * normalPDF(d1) * d_d
                                 + exp( mu_a+ pow(sd_a,2)) / sd_a  * ( -d1 ) * normalPDF(d1) * d_d
                                 - (-d2 * d_d * K * normalPDF(d2) * sd_a * S0 - K * normalPDF(d2) * sd_a )/pow(sd_a*S0,2));

        return pair<double,double>(gamma, 0.0);
    };
    pair<double,double> calcVega() const override {
        double T = m_Option.getT();
        double r = m_r;
        double v = m_Sigma;
        double t = 0;
        double S0 = m_S0;
        double K = m_K;
        double N = m_PathSize;

        double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;

        double d_mu_a = -T*((N+1.0)/(2.0*N))*(v);
        double d_sd_a = T*((N+1.0)*(2.0*N+1)/(6.0*N*N))*v/sd_a;
        double d_d2 =  (d_mu_a - d2*d_sd_a)/sd_a;
        double d_d1 = d_d2 + d_sd_a;
        double f = S0*exp(mu_a+var_a/2);
        double n1 = normalCDF(d1);
        double d_n1 = normalPDF(d1);
        double d_n2 = normalPDF(d2);

        double vega = discount(f*(n1*(d_mu_a+sd_a*d_sd_a)+ d_n1*d_d1) - K*d_n2*d_d2);

        return pair<double,double>(vega, 0.0);
    };
};

#endif //SIMULATIONMETHODS_BSASIAN_H
