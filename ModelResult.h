
#ifndef SIMULATIONMETHODS_MODELRESULT_H
#define SIMULATIONMETHODS_MODELRESULT_H

#include <ostream>
#include <string>

using namespace std;

enum class ModelType { MonteCarlo, Analytical };
enum class SensitivityMethod { FiniteDifference, PathwiseDifferentiation, LikelihoodRatio, Analytical };

class ModelResult {
    friend class Model;
    friend class MCModel;
    friend class BSModel;
    friend ostream& operator<<(ostream& os, const ModelResult &modelResult) {
        os << "Price: " << (modelResult.getPrice()) << " / Variance: " << (modelResult.getPriceVariance()) << endl;
        os << "Delta: " << (modelResult.getDelta()) << endl;
        os << "Gamma: " << (modelResult.getGamma()) << endl;
        os << "Vega: "  << (modelResult.getVega()) << endl;
        os << "Calc Time: " << (modelResult.getCalcTime()) << " ms" << endl;
        return os;
    };
    double price;
    double delta;
    double gamma;
    double vega;
    double price_variance;
    double calcTime;
    SensitivityMethod deltaMethod;
    SensitivityMethod gammaMethod;
    SensitivityMethod vegaMethod;
    ModelType modelType;
    bool antithetic;
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
    void setAntitheticVariate(bool antithetic) {
        this->antithetic = antithetic;
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
    bool usesAntitheticVariates() const { return antithetic; }
};

#endif //SIMULATIONMETHODS_MODELRESULT_H
