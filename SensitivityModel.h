
#ifndef SIMULATIONMETHODS_SENSITIVITYMODEL_H
#define SIMULATIONMETHODS_SENSITIVITYMODEL_H

#include "Path.h"
#include "Model.h"
#include "Vanilla.h"
#include "Asian.h"


/**
 * Template class responsible for helping the greeks calculation during the Monte Carlo simulations.
 * The sensitivity is calculated per path. The MCModel will then take expected value of the results from this class.
 * @tparam OPTION
 */
template<class OPTION = Option>
class SensitivityModel {
protected:
    const Model<OPTION> &m_Model;
    SensitivityMethod m_SensitivityMethod;
public:
    SensitivityModel(const Model<OPTION> &model, SensitivityMethod sensitivityMethod): m_Model(model), m_SensitivityMethod(sensitivityMethod) {}
    SensitivityMethod getSensitivityMethod() const { return m_SensitivityMethod; }
    virtual double calcDelta(const Path &path) const = 0;
    virtual double calcGamma(const Path &path) const = 0;
    virtual double calcVega(const Path &path) const = 0;
    virtual double deltaDivisor() const { return 1.0; }
    virtual double gammaDivisor() const { return 1.0; }
    virtual double vegaDivisor() const { return 1.0; }
};

namespace Greeks_by_FD {
    template<class OPTION>
    class CentralDifferencesSensitivityModel: public SensitivityModel<OPTION> {
        using SensitivityModel<OPTION>::m_Model;
    protected:
        double m_h;
    public:
        CentralDifferencesSensitivityModel(const Model<OPTION> &model, double h): SensitivityModel<OPTION>(model, SensitivityMethod::FiniteDifference), m_h(h) {}

        double getH() const { return m_h; }

        double calcDelta(const Path &path) const override {
            double payoff_bump_up = m_Model.getOption().payoff(path, Price_Up);
            double payoff_bump_down = m_Model.getOption().payoff(path, Price_Down);
            double S0 = m_Model.getS0();

            return m_h > 0 ? m_Model.discount(payoff_bump_up - payoff_bump_down): NAN;
        }

        double calcGamma(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path, None);
            double payoff_bump_up = m_Model.getOption().payoff(path, Price_Up);
            double payoff_bump_down = m_Model.getOption().payoff(path, Price_Down);
            double S0 = m_Model.getS0();

            return m_h>0 ? m_Model.discount(payoff_bump_up + payoff_bump_down - 2.0 * payoff): NAN;
        }

        double calcVega(const Path &path) const override {
            double payoff_bump_up = m_Model.getOption().payoff(path, Sigma_Up);
            double payoff_bump_down = m_Model.getOption().payoff(path, Sigma_Down);
            double S0 = m_Model.getS0();
            double sigma = m_Model.getSigma();

            return m_h>0 ? m_Model.discount(payoff_bump_up - payoff_bump_down) : NAN;
        };

        double deltaDivisor() const override { return 2 * m_h * m_Model.getS0(); }
        double gammaDivisor() const override { return m_h * m_h * m_Model.getS0() * m_Model.getS0(); }
        double vegaDivisor() const override { return 2 * m_h * m_Model.getSigma(); }
    };
}

namespace Greeks_by_PD {
    class VanillaCallSensitivityModel: public SensitivityModel<VanillaCall> {
    public:
        VanillaCallSensitivityModel(const Model<VanillaCall> &model): SensitivityModel(model, SensitivityMethod::PathwiseDifferentiation) {}
        double calcDelta(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            return payoff>0? m_Model.discount(S_T/S_0): 0.0;
        };
        double calcGamma(const Path &path) const override {
            return NAN;
        };
        double calcVega(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double sigma2 = sigma*sigma;
            double T = m_Model.getT();
            double r = m_Model.getR();
            return payoff>0 ? m_Model.discount( (S_T/sigma)*(log(S_T/S_0) - (r+0.5*sigma2)*T)): 0.0;
        };
    };

    /**
     * Pathwise Differentiation formulas for the Discrete Geometric average Fixed Strike Asian
     */
    class AsianCallSensitivityModel: public SensitivityModel<AsianCall> {
    public:
        AsianCallSensitivityModel(const Model<AsianCall> &model): SensitivityModel(model, SensitivityMethod::PathwiseDifferentiation) {}
        double calcDelta(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path,None);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            return payoff>0 ? m_Model.discount(S_T/S_0): 0.0;
        }
        double calcGamma(const Path &path) const override {
            return NAN;
        };
        double calcVega(const Path &path) const override {
            return NAN;
        };
    };
}

namespace Greeks_by_LR {
    /**
     * Can be used for options with payoffs that depend only on the terminal price, whose SDE must have an Explicit solution
     * and thus whose paths have only the terminal value.
     */
    template<class OPTION = VanillaOption>
    class VanillaSensitivityModel: public SensitivityModel<OPTION> {
        using SensitivityModel<OPTION>::m_Model;
    public:
        VanillaSensitivityModel(Model<OPTION> &model): SensitivityModel<OPTION>(model, SensitivityMethod::LikelihoodRatio) {}
        double calcDelta(const Path &path) const override {
            if(path.size()>1) return NAN; //We don't support non-Explicit SDE (path size > 1)
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double Z = path.front_random_number();
            return m_Model.discount(payoff) * Z / (S_0 * sigma * sqrt(T));
        };
        double calcGamma(const Path &path) const override {
            if(path.size()>1) return NAN; //We don't support non-Explicit SDE (path size > 1)
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double Z = path.front_random_number();
            return m_Model.discount(payoff) * ((Z*Z-1)/(S_0*S_0*sigma*sigma*T) - Z / (S_0* S_0* sigma* sqrt(T)));
        };
        double calcVega(const Path &path) const override {
            if(path.size()>1) return NAN; //We don't support non-Explicit SDE (path size > 1)
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double Z = path.front_random_number();
            return m_Model.discount(payoff) * ((Z*Z-1.0)/sigma - Z * sqrt(T));
        };
    };

    template<class OPTION = AsianOption>
    class AsianSensitivityModel: public SensitivityModel<OPTION> {
        using SensitivityModel<OPTION>::m_Model;
    public:
        AsianSensitivityModel(Model<OPTION> &model): SensitivityModel<OPTION>(model, SensitivityMethod::LikelihoodRatio) {}
        double calcDelta(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path);
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double dt = T / path.size();
            double Z = path.front_random_number();
            return m_Model.discount(payoff) * Z / sqrt(dt);
        };
        double deltaDivisor() const override { return m_Model.getS0() * m_Model.getSigma(); }
        double calcGamma(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path);
            double A_T = path.geometric_average(None);
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double score = 0.0;

            double r = m_Model.getR();
            double v = m_Model.getSigma();
            double N = path.size();

            double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
            double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
            double sd_a = sqrt(var_a);

            double Z = (log(A_T/S_0)-mu_a)/sd_a;
            return m_Model.discount(payoff) * ((Z*Z-1)/(S_0*S_0*var_a*T) - Z / (S_0* S_0* sd_a* sqrt(T)));
        };

        double calcVega(const Path &path) const override {
            double payoff = m_Model.getOption().payoff(path);
            double S_T = path.back();
            double S_0 = m_Model.getS0();
            double sigma = m_Model.getSigma();
            double T = m_Model.getT();
            double dt = T / path.size();
            double score = 0.0;
            for(int i = 0; i<path.size(); i++) {
                double Z = path.random_number(i);
                score += ((Z*Z-1.0)/sigma - Z * sqrt(dt));
            }
            return m_Model.discount(payoff) * score;
        };
    };
}

#endif //SIMULATIONMETHODS_SENSITIVITYMODEL_H
