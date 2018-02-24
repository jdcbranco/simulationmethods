//============================================================================
// Name        : Stimulation_GT.cpp
// Author      : Stanley
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <random>
using namespace std;


//model independent string comparision
bool icompare_pred(unsigned char a, unsigned char b){
		return tolower(a) == tolower(b);
}
bool icompare(string const& a, string const& b){

		if (a.length()==b.length()) {
				return equal(b.begin(), b.end(),
													 a.begin(),icompare_pred);
		}
		else {
				return false;
		}
}

//normal random variable generator
vector<double> normal_generator(unsigned int n){
		double v1;
		double v2;
		double w;
		vector<double> v;
		for(unsigned int i=0;i<(n/2)+1;i++){
				v1 = 2.0*rand()/RAND_MAX -1;
				v2 = 2.0*rand()/RAND_MAX -1;
				w = pow(v1,2) +pow(v2,2);
				if(w<=1){
					v.push_back(sqrt(-2*log(w)/w)*v1);
					v.push_back(sqrt(-2*log(w)/w)*v2);
				}else{
					--i;
				}
			}
		//size adjustment
		if(v.size()>n+1){
			v.pop_back();
			v.pop_back();
		}else{
			v.pop_back();
		}

		return v;
	}

vector<double> generate_normal(int n = 1, double mean = 0.0, double var = 1.0) {
	int j;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	normal_distribution<double> distribution(mean, var);
	vector<double> rand_numbers;
	for (j = 0; j < n; j = j + 1) {
		double random_number = distribution(generator);
		rand_numbers.push_back(random_number);
	}
	return rand_numbers;
}


void print(double d) {
	// basic print for double
	cout << d << endl;
}

void print(vector<double> &vec) {
	// basic print for vector
	for (unsigned int i = 0; i < vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

void res_print(vector<double> res) {
	cout << "price    = " << res[0] << endl;
	cout << "time     = " << res[1] << endl;
	cout << "mean     = " << res[2] << endl;
	cout << "variance = " << res[3] << endl;

}

double normalCDF(double value) {
	return 0.5 * erfc(-value / sqrt(2));
}

double normalPDF(double value) {
	return (1 / sqrt(2 * M_PI)) * exp(-0.5 * pow(value, 2));
}


class asian_option_geometric{
	  int T; // terminal time
	  unsigned int N; // number of time partitions
	  double K; // strik
	  bool is_call;


		// this function calculate pay off
		double pay_off(double g_ave){
			if ((is_call&&g_ave>K)) {
				return (g_ave - K);
			}
			else if (!is_call&&g_ave < K) {
				return -(g_ave - K);
			}
			return 0;
		}

		void MC_euler_pricing_non_anti(double S0, double r, double v, unsigned int no_sims, vector<double>& res) {
			// prameter initlisation
			clock_t t;
			double price = 0;
			double dt = (double)T / N;
			double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
			double mean_sqr = 0;

			// simulate results
			double duration;
			t = clock();// start the timer

			// run the simulation, NO variance reduction
			for (unsigned int i = 0; i < no_sims; i++) {
				vector<double> z = normal_generator(N); //generate normal vector of size N
				double log_sum = log(S0);
				double s = S0;

				for (unsigned int j = 0;j < N;j++) {
					s += r * s*dt + v * s*dt_sqr*z[j];
					log_sum += log(s);
				}
				price += pay_off(exp(log_sum / (N + 1)));
				mean_sqr += pow(pay_off(exp(log_sum / (N + 1))),2);
				
			}
			// record time
			duration = (clock() - t) / (double)CLOCKS_PER_SEC;
			// return results
			// handle edge cases
			res.push_back((price / no_sims)*exp(-r * T)); // price
			res.push_back(duration); //time
			res.push_back(price / no_sims); //mean
			res.push_back((mean_sqr / no_sims) - pow((price / no_sims), 2)); //variance

		}

	void MC_euler_pricing(double S0, double r,double v,unsigned int no_sims,vector<double>& res){
	    // prameter initlisation
	    clock_t t;
	    double price = 0;
	    double dt = (double)T/N;
	    double dt_sqr = pow((double)T/N,0.5); // ADDED: Casting T to double
		double mean_sqr = 0;

	    // simulate results
		double duration;
	    t = clock();// start the timer

	    // run the simulation and use antithetic variance reduction
	   for(unsigned int i = 0;i < no_sims;i++) {
		vector<double> z = normal_generator(N); //generate normal vector of size N
	      double sum_u = log(S0);
	      double sum_d = log(S0);
	      double s_u = S0;
	      double s_d = S0;

	      for(unsigned int j=0;j<N;j++){
			 s_u += r*s_u*dt + v*s_u*dt_sqr*z[j];
	         s_d += r*s_u*dt + v*s_u*dt_sqr*-z[j];
	         sum_u += log(s_u);
	         sum_d += log(s_d);
	      }
			// update price
			price += pay_off(exp(sum_u/(N+1)));
			price += pay_off(exp(sum_d/(N + 1)));	
			mean_sqr += pow(pay_off(exp(sum_u / (N + 1))), 2);
			mean_sqr += pow(pay_off(exp(sum_d/ (N + 1))), 2);

	    }
	    // record time
	    duration = (clock()-t)/(double)CLOCKS_PER_SEC;
		// return results
		// handle edge cases
		res.push_back((price/no_sims/2)*exp(-r*T)); // price
		res.push_back(duration); //time
		res.push_back(price/no_sims/2); //mean
		res.push_back((mean_sqr/no_sims/2)-pow((price/no_sims/2),2)); //variance

  }

	void MC_milstein_pricing(int no_sims){}
	
	void analytic_solution_pricing(double S0, double r, double v, vector<double>& res){
		clock_t c;
		double t = 0;
		// simulate results
		double duration;
		c = clock();// start the timer

		double G_t = S0;
		double mu_bar = (r - v * v / 2) * pow(T - t, 2) / (2 * T);
		double sigma_bar = sqrt(v*v / (T*T) * pow(T - t, 3) / 3);
		double d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(S0) + mu_bar - log(K));
		double d1 = d2 + sigma_bar;

		double price = exp(-r * (T - t)) * (pow(G_t, t / T) * pow(S0, (T - t) / T) * exp(mu_bar + pow(sigma_bar, 2) / 2) * normalCDF(d1) - K*normalCDF(d2));

		// record time
		duration = (clock() - c) / (double)CLOCKS_PER_SEC;

		res.push_back(price);
		res.push_back(duration);
		res.push_back(price);
		res.push_back(0.0);
	}


public :
  // an asian option has time to maturity T,K,initial price and interest rate
  asian_option_geometric(const string& type,int T,int N, double K):T(T),N(N),K(K){
    if(icompare(type,"call")){
      is_call = true;
    }
    else if(icompare(type,"put")){
      is_call = false;
    }else{
      // edge case: given comtract type known
    }
  }

  // this function calculate the option price at time 0 and analytic statistics
  vector<double> calculate_price(const string& method, double S0, double r, double v, int no_sims = 100000) {
	  // return containeer
	  vector<double> res;

	  //Use selected method to calcuate c_0
	  if (icompare(method, "euler")) {
		  MC_euler_pricing(S0, r, v, no_sims, res);
	  }
	  else if (icompare(method, "euler_na")) {
		  MC_euler_pricing_non_anti(S0, r, v, no_sims, res);
	}else if(icompare(method,"milstein")){

    }else if(icompare(method,"analytic")){
		analytic_solution_pricing(S0, r, v, res);
    }else{
      // edge case: prcing method known
    }

    return res;
  }

	// this function calculate delta at time 0 and return analytic statistics
	vector<double> calculate_delta(string method,int no_sims = 100000);

	// this function calculate gamma at time 0 and return analytic statistics
	vector<double> calculate_gamma(string method,int no_sims = 100000);


};

int main() {
	int T = 1;
	unsigned int N = 100;
	double K = 100;
	string method = "euler";
	string type = "call";
	double s0 = 100;
	double r = 0.05;
	double v = 0.4;
	int no_sims = 1000;

	srand(time(NULL));
	asian_option_geometric opt(type,T,N,K);
	vector<double> res = opt.calculate_price("analytic", s0, r, v);

	cout << "Analytic: " << endl;
	res_print(res);
	cout << endl;
	
	cout << "n = " << no_sims << endl;
	cout << "Euler: " << endl;
	res = opt.calculate_price("euler", s0, r, v, no_sims);
	res_print(res);
	cout << "Euler NA: " << endl;
	res = opt.calculate_price("euler_na", s0, r, v, no_sims);
	res_print(res);
	cout << endl;

	no_sims = 10000;
	cout << "n = " << no_sims << endl;
	cout << "Euler: " << endl;
	res = opt.calculate_price("euler", s0, r, v, no_sims);
	res_print(res);
	cout << "Euler NA: " << endl;
	res = opt.calculate_price("euler_na", s0, r, v, no_sims);
	res_print(res);
	cout << endl;

	no_sims = 100000;
	cout << "n = " << no_sims << endl;
	cout << "Euler: " << endl;
	res = opt.calculate_price("euler", s0, r, v, no_sims);
	res_print(res);
	cout << "Euler NA: " << endl;
	res = opt.calculate_price("euler_na", s0, r, v, no_sims);
	res_print(res);
	cout << endl;

	int dummy;
	cin >> dummy;
	return 0;
}
