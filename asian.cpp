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


class asian_option_geometric{
	  int T; // terminal time
	  unsigned int N; // number of time partitions
	  double K; // strik
	  bool is_call;


		// this function calculate pay off
		double pay_off(double g_ave){
		if((is_call&&g_ave>K) ){
			return (g_ave-k);
		}else if (!is_call&&g_ave < K){
			return -(g_ave-k);
		}
		return 0;
		}


	void MC_euler_pricing(double S0, double r,double v,unsigned int no_sims,vector<double>& res){
	    // prameter initlisation
	    clock_t t;
	    double price = 0;
	    double dt = (double)T/N;
	    double dt_sqr = pow(T/N,0.5);
		double mean_sqr = 0;

	    // simulate results
	    t = clock();// start the timer

	    // run the simulation and use antithetic variance reduction
	    for(unsigned int i = 0;i < no_sims/2;i++){
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
			price += pay_off(exp(sum_u*(1/(N+1))));
			price += pay_off(exp(sum_d*(1 / (N + 1))));	
			mean_sqr += pow(pay_off(exp(sum_u*(1 / (N + 1)))), 2);
			mean_sqr += pow(pay_off(exp(sum_d*(1 / (N + 1)))), 2);

	    }
	    // record time
	    t = (clock()-t)/(double)CLOCKS_PER_SEC;
		// return results
		// handle edge cases
		res.push_back((price/no_sims)*exp(-r*T)); // price
		res.push_back(t); //time
		res.push_back(price/no_sims); //mean
		res.push_back((mean_sqr/no_sims)-pow((price/no_sims),2)); //variance

  }

	void MC_milstein_pricing(int no_sims){}
	void analytic_solution_pricing(){}


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
	vector<double> calculate_price(const string& method,double S0, double r,double v,int no_sims = 100000){
    // return containeer
    vector<double> res;

    //Use selected method to calcuate c_0
    if(icompare(method,"euler")){

       	MC_euler_pricing(S0,r,v,no_sims,res);
    }else if(icompare(method,"milstein")){

    }else if(icompare(method,"analytic")){

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
	int T = 5;
	unsigned int N = 30;
	double K = 5;
	string method = "euler";
	string type = "call";
	double s0 = 10;
	double r = 0.5;
	double v = 2.3;
	int no_sims = 2;

	srand(time(NULL));
	asian_option_geometric opt(type,T,N,K);
	vector<double> res = opt.calculate_price("euler",s0,r,v,no_sims);



	return 0;
}
