
#ifndef SIMULATIONMETHODS_NORMAL_H
#define SIMULATIONMETHODS_NORMAL_H

#include <random>
#include <chrono>
#include "Random.h"

using namespace std;

class Normal: public Random {
private:
    default_random_engine generator;
protected:
    double m_Mean = 0.0;
    double m_Variance = 1.0;
    vector<double> custom_generate(unsigned int n);
    vector<double> standard_generate(unsigned int n);
public:
    Normal(const GeneratorType generatorType, double mean=0.0, double variance=1.0): Random(generatorType), generator(std::chrono::system_clock::now().time_since_epoch().count()) {
        this->m_Mean = mean;
        this->m_Variance = variance;
    }
    vector<double> generate(unsigned int n);
};


#endif //SIMULATIONMETHODS_NORMAL_H
