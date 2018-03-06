
#ifndef SIMULATIONMETHODS_MODELRESULT_H
#define SIMULATIONMETHODS_MODELRESULT_H

class ModelResult {
    friend class Model;
    friend class MCModel;
    double price;
    double delta;
    double gamma;
    double vega;
    double calcTime;
protected:
    void setPrice(double price) {
        this->price = price;
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
    double getPrice() { return this->price; }
    double getDelta() { return this->delta; }
    double getGamma() { return this->gamma; }
    double getVega() { return this->vega; }
    double getCalcTime() { return this->calcTime; }
};


#endif //SIMULATIONMETHODS_MODELRESULT_H
