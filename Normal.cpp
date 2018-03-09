#include "Normal.h"
#include <stdlib.h>
#include <chrono>
#include <random>

//This function defines which implementation to use - a custom one or a c++ standard based one
vector<double> Normal::generate(unsigned int n) {
    if(m_GeneratorType==Custom) {
        return custom_generate(n);
    } else {
        return standard_generate(n);
    }
}

//Implementation based on c++ random library, contributed by Ghali and Faycal
vector<double> Normal::standard_generate(unsigned int n) {
    int j;
    normal_distribution<double> distribution(m_Mean, m_Variance);
    vector<double> rand_numbers;
    for (j = 0; j < n; j = j + 1) {
        double random_number = distribution(generator);
        rand_numbers.push_back(random_number);
    }
    return move(rand_numbers);
}

//Implementation based on custom algorithm based on lecture notes, contributed by Stanley and Edgard
vector<double> Normal::custom_generate(unsigned int n){
    double v1;
    double v2;
    double w;
    vector<double> v;
    for(unsigned int i=0;i<(n/2)+1;i++){
        v1 = 2.0*rand()/RAND_MAX -1;
        v2 = 2.0*rand()/RAND_MAX -1;
        w = pow(v1,2)+pow(v2,2);
        if(w<=1){
            v.push_back(sqrt(-2*log(w)/w)*v1);
            v.push_back(sqrt(-2*log(w)/w)*v2);
        } else {
            --i;
        }
    }
    //size adjustment
    if(v.size()>n+1) {
        v.pop_back();
        v.pop_back();
    } else {
        v.pop_back();
    }

    return move(v);
}