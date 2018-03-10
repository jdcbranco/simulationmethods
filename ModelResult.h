
#ifndef SIMULATIONMETHODS_MODELRESULT_H
#define SIMULATIONMETHODS_MODELRESULT_H

#include <ostream>
#include <string>

using namespace std;

enum class ModelType { MonteCarlo, Analytical };
enum class SensitivityMethod { FiniteDifference, PathwiseDifferentiation, LikelihoodRatio, Analytical };

class ModelResult {
    template<typename T> friend class Model;
    template<typename T> friend class MCModel;
    template<typename T> friend class BSModel;
    template<typename T> friend class BSAsianModel;
    friend ostream& operator<<(ostream& os, const ModelResult &modelResult) {
        switch(modelResult.modelType) {
            case ModelType::MonteCarlo:
                os << "MC ";
                break;
            case ModelType::Analytical:
                os << "AN ";
                break;
        }
        switch(modelResult.deltaMethod) {
            case SensitivityMethod::Analytical:
                os << "AN ";
                break;
            case SensitivityMethod::FiniteDifference:
                os << "FD ";
                break;
            case SensitivityMethod::PathwiseDifferentiation:
                os << "PW ";
                break;
            case SensitivityMethod::LikelihoodRatio:
                os << "LR ";
        }
        os << (modelResult.usesAntitheticVariates()? "ANTI:TRUE " : "ANTI:FALSE ");
        os << (modelResult.usesControlVariate()? "CV:TRUE " : "CV:FALSE ");
        os << modelResult.getPriceVariance() << " ";
        os << sqrt(modelResult.getPriceVariance()) << " ";
        os << modelResult.getPrice() << " " << modelResult.getDelta() << " " << modelResult.getGamma() << " " << modelResult.getVega() << " " << modelResult.getCalcTime() << endl;
        return os;
    };
    double price;
    double delta;
    double gamma;
    double vega;
    double price_variance;
    double calcTime;
    ModelType modelType;
    SensitivityMethod deltaMethod;
    SensitivityMethod gammaMethod;
    SensitivityMethod vegaMethod;
    int simulations;
    bool antithetic;
    bool control_variate;
protected:
    void setPrice(double price) {
        this->price = price;
    }
    void setPriceVariance(double priceVariance) {
        this->price_variance = priceVariance;
    }
    void setDelta(double delta) {
        this->delta = delta;
    }
    void setGamma(double gamma) {
        this->gamma = gamma;
    }
    void setVega(double vega) {
        this->vega = vega;
    }
    void setCalcTime(double calcTime) {
        this->calcTime = calcTime;
    }
    void setDeltaMethod(SensitivityMethod deltaMethod) {
        this->deltaMethod = deltaMethod;
    }
    void setGammaMethod(SensitivityMethod gammaMethod) {
        this->gammaMethod = gammaMethod;
    }
    void setVegaMethod(SensitivityMethod vegaMethod) {
        this->vegaMethod = vegaMethod;
    }
    void setModelType(ModelType modelType) {
        this->modelType = modelType;
    }
    void setSimulations(int simulations) {
        this->simulations = simulations;
    }
    void setAntitheticVariate(bool antithetic) {
        this->antithetic = antithetic;
    }
    void setControlVariate(bool control_variate) {
        this->control_variate = control_variate;
    }
public:
    ModelResult() {}
    double getPrice() const { return this->price; }
    double getDelta() const { return this->delta; }
    double getGamma() const { return this->gamma; }
    double getVega() const { return this->vega; }
    double getPriceVariance() const { return this->price_variance; }
    double getCalcTime() const { return this->calcTime; }
    SensitivityMethod getDeltaMethod() const { return deltaMethod; }
    SensitivityMethod getGammaMethod() const { return gammaMethod; }
    SensitivityMethod getVegaMethod() const { return vegaMethod; }
    ModelType getModelType() const { return modelType; }
    int getSimulations() const { return this->simulations; }
    bool usesAntitheticVariates() const { return antithetic; }
    bool usesControlVariate() const { return control_variate; }
};

#endif //SIMULATIONMETHODS_MODELRESULT_H
