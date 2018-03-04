
#ifndef SIMULATIONMETHODS_NORMAL_H
#define SIMULATIONMETHODS_NORMAL_H

#include "Random.h"

class Normal: public Random {
protected:
    double m_Mean = 0.0;
    double m_Variance = 1.0;
    vector<double> custom_generate(unsigned int n);
    vector<double> standard_generate(unsigned int n);
public:
    Normal(const GeneratorType generatorType, double mean=0.0, double variance=1.0): Random(generatorType) {
        this->m_Mean = mean;
        this->m_Variance = variance;
    }
    vector<double> generate(unsigned int n);
};


#endif //SIMULATIONMETHODS_NORMAL_H
