
#ifndef SIMULATIONMETHODS_RANDOM_H
#define SIMULATIONMETHODS_RANDOM_H
#include <vector>

using namespace std;

enum GeneratorType { Custom, Standard };

class Random {
protected:
    GeneratorType m_GeneratorType;
public:
    Random(const GeneratorType generatorType): m_GeneratorType(generatorType) {}
    virtual vector<double> generate(unsigned int n) = 0;
};


#endif //SIMULATIONMETHODS_RANDOM_H
