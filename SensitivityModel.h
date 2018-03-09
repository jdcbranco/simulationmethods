//
// Created by jdcbr on 3/9/2018.
//

#ifndef SIMULATIONMETHODS_SENSITIVITYMODEL_H
#define SIMULATIONMETHODS_SENSITIVITYMODEL_H


#include "Path.h"
#include "Model.h"
#include "Vanilla.h"

class SensitivityModel {
protected:
    Model &m_Model;
    SensitivityMethod m_SensitivityMethod;
public:
    SensitivityModel(Model &model, SensitivityMethod sensitivityMethod): m_Model(model), m_SensitivityMethod(sensitivityMethod) {}
    SensitivityMethod getSensitivityMethod() const { return m_SensitivityMethod; }
    virtual double calcDelta(const Path &path) const = 0;
    virtual double calcGamma(const Path &path) const = 0;
    virtual double calcVega(const Path &path) const = 0;
};

namespace FiniteDifferences {
    class CentralDifferencesSensitivityModel: public SensitivityModel {
    public:
        CentralDifferencesSensitivityModel(Model &model): SensitivityModel(model, SensitivityMethod::FiniteDifference) {
            
        }
    };
}

namespace PathwiseDifferentiation {
    class VanillaCallSensitivityModel: public SensitivityModel {
    public:
        VanillaCallSensitivityModel(Model &model): SensitivityModel(model, SensitivityMethod::PathwiseDifferentiation) {
        }
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
}

namespace LikelihoodRatio {
    /**
     * Can be used for options with payoffs that depend only on the terminal price, whose SDE must have an Explicit solution
     * and thus whose paths have only the terminal value.
     */
    class VanillaSensitivityModel: public SensitivityModel {
    public:
        VanillaSensitivityModel(Model &model): SensitivityModel(model, SensitivityMethod::LikelihoodRatio) {
        }
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
            return m_Model.discount(payoff) * (Z*Z-1.0)/sigma - Z * sqrt(T);
        };
    };
}



#endif //SIMULATIONMETHODS_SENSITIVITYMODEL_H
