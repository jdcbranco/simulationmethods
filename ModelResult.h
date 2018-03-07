
#ifndef SIMULATIONMETHODS_MODELRESULT_H
#define SIMULATIONMETHODS_MODELRESULT_H

#include <ostream>
#include <string>

using namespace std;

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
public:
    ModelResult() {}
    double getPrice() const { return this->price; }
    double getDelta() const { return this->delta; }
    double getGamma() const { return this->gamma; }
    double getVega() const { return this->vega; }
    double getPriceVariance() const { return this->price_variance; }
    double getCalcTime() const { return this->calcTime; }

};

#endif //SIMULATIONMETHODS_MODELRESULT_H
