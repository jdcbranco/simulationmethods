#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include "Normal.h"
#include "Vanilla.h"
#include "MCModel.h"
#include "BSVanilla.h"
#include "Asian.h"

using namespace std;

Normal normal(Standard);

class Derivatives {
	double K;
	double S0;
	double T;
	double sigma;

public:
	Derivatives(double strike = 100.0, double initial_price = 100.0, double T = 1.0, double sigma = 0.4) {
		K = strike;
		S0 = initial_price;
		this->T = T;
		this->sigma = sigma;
	}

	double get_K() {
		return K;
	}
	double get_S0() {
		return S0;
	}

	double get_T() {
		return T;
	}

	double get_sigma() {
		return sigma;
	}

	void set_S0(double new_S0) {
		S0 = new_S0;
	}

	void set_sigma(double new_sigma) {
		sigma = new_sigma;
	}
};

class model{

protected:
	Derivatives d;
	double r;

public:
	model(Derivatives d, double r) {
		this->d = d;
		this->r = r;
	}

	// generate n random numbers


	double mean_vector(vector<double> vect) {
		double sum = 0.0;
		for (unsigned int i = 0; i<vect.size(); i = i + 1) {
			sum = sum + vect[i];
		}
		return sum / vect.size();
	}

	double variance_vector(vector<double> v){
		double mean = mean_vector(v);
		double sum = 0.0;
		for (int i = 0; i<v.size(); i++){
			sum += pow(v[i] - mean, 2);
		}
		sum /= (v.size() - 1);
		return pow(sum, 0.5);
	}

	virtual double CalculPrice() const {
		return 0.0;
	}

	virtual double CalculDelta() const {
		return 0.0;

	}
	virtual double CalculGamma() const {
		return 0.0;
	}

	virtual double CalculVega() const {
		return 0.0;
	}

};

class MonteCarlo : public model{
	int N;
	vector<double> price_path;
	vector<double> normal_random;

public:

	MonteCarlo(Derivatives d, double r, int N, vector<double> normal) : model(d, r) {
		this->N = N;
		this->normal_random = normal;
	}

	// generate the price path

	/*
	void euler_path() {
	int i;
	price_path = {d.get_S0()} ;
	vector<double> normal_random = normal.generate(N) ;
	for (i=1 ; i <N ; i = i+1) {
	double S_i = price_path[i-1] + r*price_path[i-1]*d.get_T()/N + d.get_sigma()*price_path[i-1]*pow(d.get_T()/N , 0.5)*normal_random[i];
	price_path.push_back(S_i);
	}
	}
	*/


	vector<double>  payoff_call(Derivatives d, int epsilon = 1) {  // Epsilon will be equal to +1 or -1 : this is to use the antithetic method for variance
		vector<double> payoff;
		int i;
		//vector<double> normal_random = normal;
		for (i = 0; i < normal_random.size(); i = i + 1) {
			double S_T = d.get_S0()*exp((r - pow(d.get_sigma(), 2.0) / 2)*d.get_T() + epsilon*d.get_sigma()*pow(d.get_T(), 0.5)*normal_random[i]);
			if (S_T - d.get_K() <= 0) {
				payoff.push_back(0.0);
			}
			else {
				payoff.push_back(exp(-r*d.get_T())*(S_T - d.get_K()));
			}
		}
		return payoff;
	}


	double variance_2(Derivatives d){
		vector<double> payoffs = payoff_call(d);
		for(int i=0;i<payoffs.size();i++){
			payoffs[i]=pow(payoffs[i],2);
		}
		double variance =mean_vector(payoffs)-(CalculPrice(d)*CalculPrice(d));
		return pow(variance,0.5);
	}

	double CalculVariance(Derivatives d){
		vector<double> payoffs = payoff_call(d);
		double variance = variance_vector(payoffs);
		return variance;

	}


	double CalculPrice(Derivatives d, int epsilon = 1) {
		vector<double> payoffs = payoff_call(d, epsilon);
		double mean_payoffs = mean_vector(payoffs);
		return mean_payoffs;
	}


	double CalculDelta(Derivatives d, double h) {
		// same derivative with price S0-h
		Derivatives d1;
		d1 = d;
		d1.set_S0(d.get_S0() - h);
		// same derivative with price S0+h
		Derivatives d2;
		d2 = d;
		d2.set_S0(d.get_S0() + h);

		// payoffs of derivatives with initial price S0-h and S0+h
		vector<double> payoff1 = payoff_call(d1);
		vector<double> payoff2 = payoff_call(d2);

		double mean_payoff1 = mean_vector(payoff1); // mean of payoffs of the derivative with initial price S0-h
		double mean_payoff2 = mean_vector(payoff2); // mean of payoffs of the derivative with initial price S0+h
		return ((mean_payoff2 - mean_payoff1) / (2 * h));

	}

	double CalculDelta_2(Derivatives d){

		vector<double> payoff1 = payoff_call(d);
		double K = d.get_K();
		double S0 = d.get_S0();
		double sum = 0;
		for (unsigned int i = 0; i<payoff1.size(); i++){
			if (payoff1[i]>0){
				sum += (payoff1[i] + K) / S0;
			}
		}
		return (sum*exp(-r*d.get_T())) / normal_random.size();


	}

	double CalculGamma(Derivatives d, double h) {
		// same derivative with price S0-h
		Derivatives d1;
		d1 = d;
		d1.set_S0(d.get_S0() - h);
		// same derivative with price S0+h
		Derivatives d2;
		d2 = d;
		d2.set_S0(d.get_S0() + h);

		// payoffs of derivatives with initial price S0-h, S0 and S0+h
		vector<double> payoff1 = payoff_call(d1);
		vector<double> payoff2 = payoff_call(d);
		vector<double> payoff3 = payoff_call(d2);

		// mean of payoffs
		double mean_payoff1 = mean_vector(payoff1); // initial price S0-h
		double mean_payoff2 = mean_vector(payoff2); // initial price S0
		double mean_payoff3 = mean_vector(payoff3); // initial price S0+h

		return ((mean_payoff1 + mean_payoff3 - 2 * mean_payoff2) / (h*h));

	}

	/*
	double CalculGamma_2(Derivatives d){
		vector<double> payoff1 = payoff_call(d);
		double K = d.get_K();
		double S0 = d.get_S0();
		double sum = 0;
		double T = d.get_T();
		double sigma = d.get_sigma();
		double d1 = (log(S0 / K) + ((r + pow(sigma, 2) / 2))*T)/(sigma*pow(T,0.5));
		for (unsigned int i = 0; i<payoff1.size(); i++){
			if (payoff1[i]>0){
				sum += (payoff1[i] + K) ;
			}
		}
		return (sum*exp(-r*T)*(pow(d1,2.0)-d1*sigma*pow(T,0.5)-1) /(S0*S0*sigma*sigma*T* normal_random.size()));

	}
	*/
	double CalculVega(Derivatives d, double h) {
		// same derivative with volatility sigma-h
		Derivatives d1;
		d1 = d;
		d1.set_sigma(d.get_sigma() - h);
		// same derivative with volatility sigma+h
		Derivatives d2;
		d2 = d;
		d2.set_sigma(d.get_sigma() + h);

		// payoffs of derivatives with initial volatility sigma-h and sigma+h
		vector<double> payoff1 = payoff_call(d1);
		vector<double> payoff2 = payoff_call(d2);

		double mean_payoff1 = mean_vector(payoff1); // mean of payoffs of the derivative with initial volatility sigma-h
		double mean_payoff2 = mean_vector(payoff2); // mean of payoffs of the derivative with initial volatility sigma+h
		return ((mean_payoff2 - mean_payoff1) / (2 * h));
	}


	double CalculVega_2(Derivatives d){
		vector<double> payoff1 = payoff_call(d);
		double K = d.get_K();
		double S0 = d.get_S0();
		double sum = 0;
		double sigma = d.get_sigma();
		double ST;
		double T = d.get_T();
		for (unsigned int i = 0; i<payoff1.size(); i++){
			if (payoff1[i]>0){
				ST = payoff1[i] + K;
				sum += ((ST) / sigma)*(log(ST / S0) - (r + 0.5*pow(sigma, 2)*T));
			}
		}
		return (sum*exp(-r*d.get_T())) / normal_random.size();
	}


};


class Black_Scholes :public model{

public:
	Black_Scholes(Derivatives d, double r) : model(d, r) { }

	double normalCDF(double value) {
		return 0.5 * erfc(-value / sqrt(2));
	}

	double CalculPrice() {
		double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
		double d2 = d1 - d.get_sigma() * sqrt(d.get_T());
		return (d.get_S0()*normalCDF(d1) - normalCDF(d2)*d.get_K()*exp(-r*d.get_T()));
	}

	double CalculDelta() {
		double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
		return (normalCDF(d1));
	}

	double CalculGamma() {
		double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
		return (exp(-d1*d1 / 2) / sqrt(2 * M_PI) / (d.get_S0()*d.get_sigma()*sqrt(d.get_T())));
	}

	double CalculVega() {
		double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
		return(exp(-d1*d1 / 2) / sqrt(2 * M_PI) * d.get_S0() *sqrt(d.get_T()));
	}
};


vector<double> generate(const Derivatives &d, int number_simulations , double r = 0.05) {
    vector<double> results;
    vector<double> normal_vector = normal.generate(number_simulations);
    MonteCarlo mc(d, r, 100, normal_vector);
    double price_call_mc = mc.CalculPrice(d);
    double variance_mc = mc.variance_2(d);
    double delta_mc = mc.CalculDelta(d, 0.1);
    double gamma_mc = mc.CalculGamma(d,0.5);
    double vega_mc = mc.CalculVega(d, 0.5);
    results.push_back(number_simulations);
    results.push_back(price_call_mc);
    results.push_back(variance_mc);
    results.push_back(delta_mc);
    results.push_back(gamma_mc);
    results.push_back(vega_mc);
    return results;
}

double mean_vector(vector<double> vect) {
    double sum = 0.0;
    for (unsigned int i = 0; i<vect.size(); i = i + 1) {
        sum = sum + vect[i];
    }
    return sum / vect.size();
}

double variance_vector(vector<double> v){
    double mean = mean_vector(v);
    double sum = 0.0;
    for (int i = 0; i<v.size(); i++){
        sum += pow(v[i] - mean, 2);
    }
    sum /= (v.size() - 1);
    return pow(sum, 0.5);
}

vector<double> generate_mean_std(const Derivatives &d, int total_simulations,int n_simulations,double r = 0.05) {
    // total_simulations : total number of call simulations
    // n_simulations : number of simulation to generate a call price
    int i;
    vector<double> results;
    for (i=0; i<total_simulations;i=i+1) {
        vector<double> normal_vector = normal.generate(n_simulations);
        MonteCarlo mc(d, r, 100, normal_vector);
        double price_call_mc = mc.CalculPrice(d);
        results.push_back(price_call_mc);
    }
    double mean = mean_vector(results);
    double std = variance_vector(results);
    vector<double> res = {mean,std};
    return res;
}

int main() {
	double strike = 100;
	//int number_simulations = 100000;
	double sigma = 0.4;
	double r = 0.05;

	vector<int> number_simulations = {1000,5000,10000,25000,50000,100000,200000,500000} ;

    //New code for European call
    VanillaCall vanillaCall(strike, 1.0);
    BSCallModel bsModel(vanillaCall, 100.0, sigma, r);
    MCModel mcModel(vanillaCall, 100.0, sigma, r, 0.005, Explicit);
    Simulator simulator(normal,true); //True=Antithetic
    ModelResult bsModelResult = bsModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes: " << endl;
    cout << bsModelResult;
	for(int i: number_simulations) {
		ModelResult mcModelResult = mcModel.simulate(simulator,i,mcModel.getSolver() == Explicit ? 1 : 5);
		cout << "-------------" << endl;
		cout << "Simulations: " << i << endl;
		cout << mcModelResult;
	}
    //New code for Asian call
    print_sample_path = true;
    GFixedStrikeAsianCall asianCall(strike, 1.0);
    MCModel asianMcModel(asianCall, 100.0, sigma, r, 0.005, Explicit); //Optionally, can try Euler as well. Both work fine.
    ModelResult asianMcModelResult = asianMcModel.simulate(simulator,100000, 100);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and 100 steps: " << endl;
    cout << asianMcModelResult;

	//End-New
    return 0;

	Derivatives call(strike, 100, 1.0, sigma);

    
    ofstream fout("output.txt"); // creates an ofstream called fout
    if (! fout.is_open()) { // test that file is open
        cout << "Error opening output file." << endl;
        return -1;
    }
    
    for (int i = 0 ; i<number_simulations.size(); i=i+1) {
        // fout << i+1 << "," ;
        vector<double> res_sim = generate(call, number_simulations[i]);
        fout << res_sim[0] ;
        for (int j = 1 ; j < res_sim.size() ; j=j+1) {
            fout << "," << res_sim[j]   ;
        }
        fout << endl;
    }
    fout.close();
    
    
    vector<double> total_sim = {1000} ;
    vector<double> sim_per_call = {1000,10000,20000,50000,100000,500000};
    
    
    ofstream fout2("output_2.txt"); // creates an ofstream called fout
    if (! fout2.is_open()) { // test that file is open
        cout << "Error opening output file." << endl;
        return -1;
    }
    
    for (int j = 0 ; j < total_sim.size(); j=j+1) {
        for (int k = 0; k< sim_per_call.size();k=k+1) {
            vector<double> results = generate_mean_std(call, total_sim[j], sim_per_call[k]);
            fout2 << total_sim[j] << " " << sim_per_call[k] << " " << results[0] << " " << results[1] << endl;
        }
    }

    fout2.close();

	std::clock_t    start;
	start = std::clock();

	vector<double> normal_vector = normal.generate(number_simulations[2]);
	std::cout << "Time taken to compute the normal vector: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	MonteCarlo mc(call, r, 100, normal_vector);
	Black_Scholes bs(call, r);


	// price of the option
	start = std::clock();
	double price_call_mc = mc.CalculPrice(call);
	cout << "The Monte_Carlo price of the call is equal to " << price_call_mc << endl;
	std::cout << "Time taken to compute the Monte_Carlo price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	start = std::clock();
	double price_call_bs = bs.CalculPrice();
	cout << "The BS price of the call is equal to " << price_call_bs << endl;
	std::cout << "Time taken to compute the BS price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;


	double price_call_MC_antit = (mc.CalculPrice(call) + mc.CalculPrice(call ,-1)) / 2;
	cout << "The MC price of the call using antithetic variables is equal to " << price_call_MC_antit << endl;

	vector<double> v = mc.payoff_call(call, +1);

	vector<double> call_price_vector;

	for (int i=0;i<1000;i++){
		vector<double> normal_vector = normal.generate(number_simulations[2]);
		MonteCarlo mc(call, r, 100, normal_vector);
		call_price_vector.push_back(mc.CalculPrice(call));
	}
	double variance_newmethod = mc.variance_vector(call_price_vector);
	cout << "The Monte_Carlo variance of the call using the new method is  " << variance_newmethod << endl;

	//double variance_antiT = mc.variance_vector(v);
	//cout << "The Monte_Carlo variance of the call using Antitetic variables is equal to " << variance_antiT << endl;


	//mean

	double variance_MC = mc.variance_2(call);
	cout << "The Monte_Carlo variance of the call is equal to " << variance_MC << endl;


	// delta
	/*double delta_mc = mc.CalculDelta(call, 1000000, 0.1);
	cout << "The Monte_Carlo delta is equal to "<< delta_mc << endl;*/
	start = std::clock();
	double delta_bs = bs.CalculDelta();
	std::cout << "Time taken to compute the delta_BS price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	cout << "The BS delta is equal to " << delta_bs << endl;

	start = std::clock();
	double delta_new = mc.CalculDelta_2(call);
	std::cout << "Time taken to compute the delta_1 MC price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	cout << "The Monte_Carlo delta_1 is equal to " << delta_new << endl;

	// gamma

	start = std::clock();
	double gamma_mc = mc.CalculGamma(call,0.5);
	cout << "The MonteCarlo gamma is equal to " << gamma_mc << endl;
	std::cout << "Time taken to compute the Gamma MC price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	start = std::clock();
	double gamma_bs = bs.CalculGamma();
	cout << "The BS gamma is equal to " << gamma_bs << endl;
	std::cout << "Time taken to compute the BS Gamma price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	// vega
	start = std::clock();
	double vega_mc = mc.CalculVega_2(call);
	std::cout << "Time taken to compute the MC Vega price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	cout << "The Monte Carlo vega is equal to " << vega_mc << endl;

	start = std::clock();
	double vega_bs = bs.CalculVega();
	cout << "The BS vega is equal to " << vega_bs << endl;

	std::cout << "Time taken to compute the BS Vega price: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	return 0;
}
