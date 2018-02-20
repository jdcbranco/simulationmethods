#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
#include <chrono>
using namespace std;


class Derivatives {
    double K;
    double S0;
    
public:
    Derivatives(double strike, double initial_price) {
        K = strike;
        S0 = initial_price;
    }
    
    double get_K() {
        return K;
    }
    double get_S0() {
        return S0;
    }
    
    void set_S0(double new_S0) {
        S0 = new_S0;
    }
    
};
    
    

class model{
    
protected:
    Derivatives* d;
    double r;
    double T;
    double sigma;
    
public :
    model (Derivatives* d,double r,double T,double sigma) {
        this ->d = d;
        this->r = r;
        this-> T = T;
        this -> sigma =sigma ;
    }
    
    // generate n random numbers
    vector<double> generate_normal(double mean = 0.0 ,double var = 1.0 ,int n=1) {
        int j ;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        normal_distribution<double> distribution(mean,var) ;
        vector<double> rand_numbers ;
        for (j = 0 ; j < n ; j = j+1) {
            double random_number = distribution(generator);
            rand_numbers.push_back(random_number) ;
        }
        return rand_numbers;
    }
    
    
    virtual double CalculPrice() const =0 ;
    virtual double CalculDelta() const =0;
    virtual double CalculGamma() const =0;
    
};


class MonteCarlo : public model{
    int N;
    vector<double> price_path;
    
public :
    
    MonteCarlo (Derivatives* d,double r,double T,double sigma, int N) : model(d,r,T,sigma) {
        this -> N = N;
    }
    
    // generate the price path
    void euler_path() {
        int i;
        price_path = {d -> get_S0()} ;
        vector<double> normal_random = generate_normal(0.0, 1.0, N) ;
        for (i=1 ; i <N ; i = i+1) {
            double S_i = price_path[i-1] + r*price_path[i-1]*T/N + sigma*price_path[i-1]*pow(T/N , 0.5)*normal_random[i];
            price_path.push_back(S_i);
        }
    }
    
    vector<double> payoff_call(Derivatives* d, int n_simulations) {
        vector<double> payoff ;
        int i ;
        vector<double> normal_random = generate_normal(0.0, 1.0, n_simulations) ;
        for (i = 0 ; i < n_simulations ; i = i+1) {
            double S_T = d -> get_S0()*exp( (r- pow(sigma,2.0)/2)*T + sigma*pow(T, 0.5)*normal_random[i]);
            if (S_T-d -> get_K() < 0 ) {
                payoff.push_back(0.0)  ;
            }
            else {
                payoff.push_back(exp(-r*T)*(S_T-d -> get_K()));
            }
        }
        return payoff;
    }
    
    double CalculDelta(Derivatives* d, int n_simulations,double h) {
        int i;
        // same derivative with price S0-h
        Derivatives* d1;
        *d1 = *d;
        d1 -> set_S0(d -> get_S0() -h) ;
        // same derivative with price S0+h
        Derivatives* d2;
        *d2 = *d;
        d2 -> set_S0(d -> get_S0() +h) ;
        
        // payoff stores payoffs of the same derivative
        vector<double> payoff1 = payoff_call(d1, n_simulations) ;
        vector<double> payoff2 = payoff_call(d2, n_simulations) ;
        
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (i =0 ; i<n_simulations; i=i+1) {
            sum1 = sum1 + payoff1[i];
            sum2 = sum2 + payoff2[i];
        }
        return ((sum2-sum1)/(2*h*n_simulations));
        
    }
    

};


class Black_Scholes :public model{
    
public:
    Black_Scholes(Derivatives* d,double r,double T,double sigma) : model(d,r,T,sigma) { }
    
    double normalCDF(double value) {
        return 0.5 * erfc(-value / sqrt(2));
    }
    
    double CalculPrice() {
        double d1 = (log(d->get_S0() / d->get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        return (d->get_S0()*normalCDF(d1) - normalCDF(d2)*d->get_K()*exp(-r*T) );
    }
    
    double CalculDelta() {
        double d1 = (log(d->get_S0() / d->get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        return (normalCDF(d1));
    }
                
    double CalculGamma() {
        double d1 = (log(d->get_S0() / d->get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        return (exp(-d1*d1/2)/sqrt(2*M_PI) /(d->get_S0()*sigma*sqrt(T)) ) ;
    }
    
};





