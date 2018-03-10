#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "Normal.h"
//#define M_PI 3.141592653589793238462643383279502884

using namespace std;
void hline(int n) {
    for(int i = 0; i < n; i++ ) {
        cout << "-";
    }
    cout << endl;
}

void dotline(int n) {
    for(int i = 0; i < n; i++ ) {
        cout << ".";
    }
    cout << endl;
}


void output_file(string file_name, vector<double> &vec ) {
    ofstream fout;
    fout.open (file_name);
    for(unsigned int i = 0; i < vec.size(); i++ ) {
        fout << vec[i] << endl;
    }
    fout.close();
}

void output_file(string file_name, vector< vector<double> > &mat ) {
    ofstream fout;
    fout.open (file_name);
    for(unsigned int i = 0; i < mat.size(); i++ ) {
        for(unsigned int j = 0; j < mat[i].size() - 1; j++ ) {
            fout << mat[i][j] << ",";
        }
        fout << mat[i][mat[i].size()-1];
        fout << endl;
    }
    fout.close();
}


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

void res_print(vector<double> &res) {

    if(res.size() < 3) {
        cout << "   value    = " << res[0] << endl;
        cout << "   time     = " << res[1] << endl;
    } else {
        cout << "   estimate = " << res[0] << endl;
        cout << "   time     = " << res[1] << endl;
        cout << "   err_mean = " << res[2] << endl;
        cout << "   err_var  = " << res[3] << endl;
    }
}

double normalCDF(double value) {
	return 0.5 * erfc(-value / sqrt(2));
}

double normalPDF(double value) {
	return (1 / sqrt(2 * M_PI)) * exp(-0.5 * pow(value, 2));
}



/******************************************************************************/

//Normal random number generator
Normal normal(Custom);

class asian_option_geometric{
    //################## why T is a double ################
    int  T; // terminal time
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

    void MC_euler_pricing(double S0, double r, double v, unsigned int no_sims, vector<double>& res) {
        // prameter initlisation
        clock_t c;
        double duration;
        double dt = (double)T / N;
        double dt_sqrt = pow(dt, 0.5); // ADDED: Casting T to double
        double this_price = 0; // itermediate price results
        double price = 0; // price
        double discount = exp(-r * T); // discouting factor
        double log_sum = 0; // intemediate log price sum
        double s = S0; // intemediate pices
        // run the simulation, NO variance reduction
        c = clock();
        for (unsigned int i = 0; i < no_sims; i++) {
            vector<double> z = normal.generate(N); //generate normal vector of size N

            for (unsigned int j = 0;j < N;j++) {
                s += r * s*dt + v * s*dt_sqrt*z[j];
                log_sum += log(s);
            }
            this_price = pay_off(exp(log_sum / N));
            price += this_price;
            log_sum = 0; // rets path
            s = S0; // rets path
        }

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        price = (price * discount/ no_sims);

        res.push_back(price); // price
        res.push_back(duration); //time
    }

		void MC_analytic_pricing(double S0, double r, double v, unsigned int no_sims, vector<double>& res){
				clock_t c;
				double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
				double sd_a = sqrt(var_a);
				double discount = exp(-r * T);
				double price = 0;
				c = clock();// start the timer

				// run simulation
				vector<double> z = normal.generate(no_sims);
				for (unsigned int i = 0; i < no_sims; i++) {
					price += pay_off((S0) * exp(mu_a + sd_a *z[i]));
				}

				// calculate detla
				price = discount*price/no_sims;

				// record time
				duration = (clock() - c) / (double)CLOCKS_PER_SEC;

				res.push_back(price);
				res.push_back(duration);
		}

		void analytic_solution_pricing(double S0, double r, double v, vector<double>& res){
        double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;
				double discount  = exp(-r*T);
        double price = 0;
				if(is_call){
					price = discount*(S0 * exp(mu_a + 0.5*var_a)*normalCDF(d1) - K*normalCDF(d2));
				}else{
				}

        // return results
        res.push_back(price);
        res.push_back(0);
    }

    void MC_fd_delta(double S0, double r, double v,unsigned int no_sims, vector<double> &res, double h) {
        clock_t c;
        double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double discount = exp(-r * T);
        double price1 = 0;
				double price2 = 0;
        double delta = 0;
        c = clock();// start the timer

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for (unsigned int i = 0; i < no_sims; i++) {
          price1 += pay_off((S0+h) * exp(mu_a + sd_a *z[i]));
          price2 += pay_off((S0-h) * exp(mu_a + sd_a *z[i]));
        }

        // calculate detla
        price1 = discount*price1/no_sims;
        price2 = discount*price2/no_sims;
				delta = (price1-price2)/(2*h);

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;

        res.push_back( delta );
        res.push_back( duration );

    }

    void MC_pw_delta(double S0, double r, double v, unsigned int no_sims, vector<double> &res) {
        clock_t c;
        double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double discount = exp(-r * T);
        double delta = 0;
        double payoff = 0; // intemediate results
        // run simulation
        c = clock();// start the timer
        vector<double> z = normal.generate(no_sims);
        for (unsigned int i = 0; i < no_sims; i++) {
            payoff = pay_off(S0 * exp(mu_a + sd_a *z[i]));
            if(payoff>0 && is_call){
                delta += exp(mu_a + sd_a *z[i]);
            } else if (payoff>0 && !is_call) {
                delta += -exp(mu_a + sd_a *z[i]);
            }
        }

        // calculate detla
        delta = discount*delta/no_sims;
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        // record results
        res.push_back( delta );
        res.push_back( duration );

    }

    void MC_lr_delta(double S0, double r, double v, unsigned int no_sims, vector<double> &res){
        clock_t c;
        double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double delta = 0;
				double discount = exp(-r*T);


        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(unsigned int i=0;i<no_sims;i++){
          delta += pay_off(S0 * exp(mu_a + sd_a*z[i]))*z[i]/(sd_a*S0);
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        delta = discount*delta/no_sims;

        // return the display_results
        res.push_back(delta); //delta
        res.push_back(duration); //time


    }

		void analytic_solution_delta(double S0, double r,double v, vector<double> &res){
				double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
				double sd_a = sqrt(var_a);
				double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
				double d1 = d2+sd_a;
				double discount  = exp(-r*T);
				double delta = 0;
				if(is_call){
							delta = discount*(exp(mu_a+0.5*var_a)*(normalCDF(d1)+normalPDF(d1)/sd_a) - K*normalPDF(d2)/(sd_a*S0));
				}else{

				}
		    res.push_back(delta);
		    res.push_back(0);
			}

    void MC_fd_gamma(double S0, double r, double v, unsigned int no_sims, vector<double> &res, double h) {

      clock_t c;
      double duration;
			double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
			double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
			double sd_a = sqrt(var_a);
      double discount = exp(-r * T);
      double price1 = 0;
			double price2 = 0;
			double price3 = 0;
      double gamma = 0;
      c = clock();// start the timer

      // run simulation
      vector<double> z = normal.generate(no_sims);
      for (unsigned int i = 0; i < no_sims; i++) {
        price1 += pay_off((S0+h) * exp(mu_a + sd_a *z[i]));
        price2 += pay_off(S0 * exp(mu_a + sd_a *z[i]));
        price3 += pay_off((S0-h) * exp(mu_a + sd_a *z[i]));
      }

      // calculate detla
      price1 = discount*price1/no_sims;
      price2 = discount*price2/no_sims;
      price3 = discount*price3/no_sims;
      gamma = (price1-2*price2+price3)/(h*h);

      // record time
      duration = (clock() - c) / (double)CLOCKS_PER_SEC;

      res.push_back(gamma);
      res.push_back(duration);
    }

    void MC_lr_gamma(double S0, double r, double v, unsigned int no_sims, vector<double> &res){
        clock_t c;
        double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
        double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
        double sd_a = sqrt(var_a);
        double S0_sqr = S0*S0;
        double gamma = 0;
				double discount = exp(-r*T);

        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(unsigned int i=0;i<no_sims;i++){
          gamma += pay_off(S0 * exp(mu_a + sd_a*z[i]))*((z[i]*z[i]-1)/(S0_sqr*var_a)-z[i]/(S0_sqr*sd_a));
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        gamma = discount*gamma/no_sims;

        // return the display_results
        res.push_back(gamma); //gamma
        res.push_back(duration); //time
    }

		void analytic_solution_gamma(double S0, double r, double v, vector<double> &res){
				double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
				double sd_a = sqrt(var_a);
				double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
				double d1 = d2+sd_a;
				double discount  = exp(-r*T);
				double d_d = 1.0 / (sd_a*S0);
				double gamma = 0;
				if(is_call){
					gamma =  discount*(exp( mu_a + pow(sd_a,2 ) ) * normalPDF(d1) * d_d
				                    + exp( mu_a+ pow(sd_a,2)) / sd_a  * ( -d1 ) * normalPDF(d1) * d_d
				                    - (-d2 * d_d * K * normalPDF(d2) * sd_a * S0 - K * normalPDF(d2) * sd_a )/pow(sd_a*S0,2));
				}else{
				}
				res.push_back(gamma);
		    res.push_back(0);
    }

		void MC_fd_vega(double S0, double r, double v, unsigned int no_sims, vector<double> &res, double h) {
      clock_t c;
      double duration;
      double vega = 0;
      // prameters for the first price
      double mu_a1 = T*(r-(v+h)*(v+h)/2)*((N+1.0)/(2.0*N));
      double var_a1 = T*(v+h)*(v+h)*(1.0/N+(N-1.0)*(2.0*N-1)/(6.0*N*N));
      double sd_a1 = sqrt(var_a1);

      // prameters for the second price
      double mu_a2 = T*(r-(v-h)*(v-h)/2)*((N+1.0)/(2.0*N));
      double var_a2 = T*(v-h)*(v-h)*(1.0/N+(N-1.0)*(2.0*N-1)/(6.0*N*N));
      double sd_a2 = sqrt(var_a2);
      double discount = exp(-r * T);
      double price1 = 0;
			double price2 = 0;

      c = clock();// start the timer

      // run simulation
      vector<double> z = normal.generate(no_sims);
      for (unsigned int i = 0; i < no_sims; i++) {
        price1 += pay_off(S0 * exp(mu_a1 + sd_a1 *z[i]));
        price2 += pay_off(S0 * exp(mu_a2 + sd_a2 *z[i]));
      }

      // calculate detla
      price1 = discount*price1/no_sims;
      price2 = discount*price2/no_sims;
      vega = (price1-price2)/(2*h);

      // record time
      duration = (clock() - c) / (double)CLOCKS_PER_SEC;

      res.push_back( vega );
      res.push_back( duration );
    }

		void MC_pw_vega(double S0, double r, double v, unsigned int no_sims, vector<double> &res){
	      clock_t c;
	      double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*(1.0/3 - 1.0/(2*N) + 1.0/(6.0*N));
				double sd_a = sqrt(var_a);
      	double d_mu_a = -v*T*(1.0/N +(N-1.0)/(2.0*N));
				double d_sd_a = T*v*(1.0/N+(N-1.0)*(2.0*N-1)/(6.0*pow(N,2)))/sd_a;
				double vega = 0;
				double A = 0; // record inte results;
				double payoff = 0;
				double discount = exp(-r*T);


	      // start the timer
	      c = clock();

	      // run simulation
	      vector<double> z = normal.generate(no_sims);
	      for(unsigned int i=0;i<no_sims;i++){
					 A = S0 * exp(mu_a + sd_a *z[i]);
					 payoff = pay_off(A);
					if((payoff>0 && is_call) || (payoff>0 && !is_call)){
						vega += A*(d_mu_a +z[i]*d_sd_a);
					}
	      }

	      // stop the timer
	      duration = (clock()-c)/(double)CLOCKS_PER_SEC;

	      // calculate the results
	      vega = vega * discount/no_sims;

	      // return the display_results
	      res.push_back(vega); //vega
	      res.push_back(duration); //time
    }

		void MC_lr_vega(double S0, double r, double v, unsigned int no_sims, vector<double> &res){
        clock_t c;
        double duration;
				double mu_a = T*(r-v*v/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*(1.0/3 - 1.0/(2*N) + 1.0/(6.0*N));
				double sd_a = sqrt(var_a);
      	double d_mu_a = -v*T*(1.0/N +(N-1.0)/(2.0*N));
				double d_sd_a = T*v*(1.0/N+(N-1.0)*(2.0*N-1)/(6.0*pow(N,2)))/sd_a;
				double g = d_sd_a/sd_a;
				double df_f = 0;

				double vega = 0;
  			double discount = exp(-r*T);

        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(unsigned int i=0;i<no_sims;i++){
					df_f = (z[i])*(d_mu_a+z[i]*d_sd_a)/sd_a - g;
					vega += pay_off(S0 * exp(mu_a + sd_a*z[i]))*df_f;
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        vega = discount*vega/no_sims;

        // return the display_results
        res.push_back(vega); //vega
        res.push_back(duration); //time
    }

		void analytic_solution_vega(double S0, double r, double v, vector<double>& res){
				double mu_a = T*(r-pow(v,2)/2)*((N+1.0)/(2.0*N));
				double var_a = T*v*v*((N+1.0)*(2.0*N+1)/(6.0*pow(N,2)));
				double sd_a = sqrt(var_a);
				double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
				double d1 = d2+sd_a;
                double d_mu_a = -T*((N+1.0)/(2.0*N))*(v);
                double d_sd_a = T*((N+1.0)*(2.0*N+1)/(6.0*N*N))*v/sd_a;
                double d_d2 =  (d_mu_a - d2*d_sd_a)/sd_a;
				double d_d1 = d_d2 + d_sd_a;
            double f = S0*exp(mu_a+var_a/2);
                double n1 = normalCDF(d1);
                double d_n1 = normalPDF(d1);
                double d_n2 = normalPDF(d2);
                double discount = exp(-r*T);
            
				double vega = 0;
				if(is_call){
                    vega =  discount*(f*(n1*(d_mu_a+sd_a*d_sd_a)+ d_n1*d_d1) - K*d_n2*d_d2);
				}else{}
				// record time
				res.push_back(vega);
				res.push_back(0);
    }

    void display_opening(string target_val, string method1, string method2 ) {
        hline(50);
        cout << target_val <<" using method: " << method1 << " with " << method2 << endl;
    }

    void display_opening(string target_val, string method1) {
        hline(50);
        cout << target_val <<" using method: " << method1 << endl;
    }

    void display_ending() {
        hline(50);
    }

    void display_section_split() {
        dotline(50);
    }
  public:
    // an asian option has time to maturity T,K,initial price and interest rate
    asian_option_geometric(const string& type,int T,int N, double K):T(T),N(N),K(K){
        if(icompare(type,"call")){
            is_call = true;
        }
        else if(icompare(type,"put")){
            is_call = false;
        } else {
            // edge case: given comtract type known
        }
    }

    // this function calculate the option price at time 0 and analytic statistics
    vector<double> calculate_price(const string& method, double S0, double r, double v, unsigned int no_sims = 100000, bool display_results = false) {
        // return containeer
        vector<double> res;

        //Use selected method to calcuate c_0
        if (icompare(method, "euler")) {
            MC_euler_pricing(S0, r, v, no_sims, res);
        } else if (icompare(method,"analytic")){
            analytic_solution_pricing(S0, r, v, res);
        } else if (icompare(method,"mc analytic")){
            MC_analytic_pricing(S0, r, v, no_sims,res);
				} else {
            cout << "method " << method << " not recognized, please choose another" << endl;
        }

        string name = "PRICE";
        if(display_results) {
            display_opening(name,method);

            display_section_split();
						cout << "   K       = " << K << endl;
						cout << "   T       = " << T << endl;
            cout << "   S0       = " << S0 << endl;
            cout << "   r        = " << r  << endl;
            cout << "   v        = " << v  << endl;

            if( !icompare(method,"analytic") ) {
								cout << "   N       = " << N << endl;
                cout << "   M        = " << no_sims << endl;
            }
            display_section_split();
            res_print(res);
            display_ending();
        }

        return res;

    }

    // this function calculate delta at time 0 and return analytic statistics
    vector<double> calculate_delta(string method2, double S0, double r, double v, unsigned int no_sims = 100000, double h = 0.01, bool display_results = false ) {
        // method is between euler/ emilistein
        // method2 is the user choice of between finite difference/ pathwise/ likelihood ratio
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_delta(S0, r, v, no_sims, res, h);

        } else if(icompare(method2, "pw")) {
            MC_pw_delta(S0, r, v, no_sims, res);


        } else if(icompare(method2, "analytic") ) {
            analytic_solution_delta(S0,r,v,res);

				}else if(icompare(method2, "lr") ) {
					  MC_lr_delta(S0,r,v,no_sims,res);
			  }else {
            // edge case:
        }

        string name = "DELTA";
        if(display_results) {
            display_opening(name,method2);

            display_section_split();

            cout << "   S0       = " << S0 << endl;
            cout << "   r        = " << r  << endl;
            cout << "   v        = " << v  << endl;

            if( !icompare(method2,"analytic") ) {
                cout << "   M        = " << no_sims << endl;
                cout << "   h        = " << h << endl;
            }

            display_section_split();

            res_print(res);
            display_ending();
        }

        return res;
    };

		//this function calculate vega at time 0 and return analytic statistics
    vector<double> calculate_vega(string method2, double S0, double r, double v, unsigned int no_sims = 100000, double h = 0.01, bool display_results = false) {
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_vega(S0, r, v, no_sims, res,h);
        } else if(icompare(method2, "analytic")){
            analytic_solution_vega(S0,r,v,res);
        }else if(icompare(method2, "lr")){
          	MC_lr_vega(S0,r,v,no_sims,res);
				} else if(icompare(method2, "pw")){
          	MC_pw_vega(S0,r,v,no_sims,res);
				}else{
            // edge case:
        }

        string name = "VEGA";
        if(display_results) {
            display_opening(name,method2);

            display_section_split();

            cout << "   S0       = " << S0 << endl;
            cout << "   r        = " << r  << endl;
            cout << "   v        = " << v  << endl;

            if( !icompare(method2,"analytic") ) {
                cout << "   M        = " << no_sims << endl;
                cout << "   h        = " << h << endl;
            }

            display_section_split();

            res_print(res);
            display_ending();
        }

        return res;
    };

    // this function calculate vega at time 0 and return analytic statistics
    vector<double> calculate_gamma(string method2, double S0, double r, double v, unsigned int no_sims = 100000, double h = 0.01, bool display_results = false) {
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_gamma(S0,r, v, no_sims, res, h);
        } else if(icompare(method2, "analytic")) {
            analytic_solution_gamma(S0,r,v,res);
        }else if(icompare(method2, "lr")) {
            MC_lr_gamma(S0,r,v,no_sims,res);
        } else {
            // edge case:
        }

        string name = "GAMMA";
        if(display_results) {
            display_opening(name,method2);

            display_section_split();

            cout << "   S0       = " << S0 << endl;
            cout << "   r        = " << r  << endl;
            cout << "   v        = " << v  << endl;

            if( !icompare(method2,"analytic") ) {
                cout << "   M        = " << no_sims << endl;
                cout << "   h        = " << h << endl;
            }

            display_section_split();

            res_print(res);
            display_ending();
        }


        return res;
    }
};


int main() {
	int T = 1;
	unsigned int N = 1000;
	double K = 100;
	string method = "euler";
	string type = "call";
	double s0 = 100;
	double r = 0.05;
	double v = 0.4;
	//unsigned int no_sims = 100000;
    double h = 0.01;
    int no_sims_lower = 1000;
    int no_sims_upper = 110000;
    int no_sims_step = 1000;
    
    // create an asian geometric average
    asian_option_geometric opt("call",T,N,K);

//    opt.calculate_vega("analytic",s0,r,v,false);
//    opt.calculate_vega("pw",s0,r,v,no_sims,false);
//    opt.calculate_vega("fd",s0,r,v,no_sims,h,false);
//    opt.calculate_vega("lr",s0,r,v,no_sims,false);
//
//    opt.calculate_delta("pw",s0,r,v,no_sims,h,true);
//    opt.calculate_delta("lr",s0,r,v,no_sims,h,true);

    
     // write files
     ofstream myfile;
     myfile.open("prices vs no_sim.csv");
     ofstream myfile1;
     myfile1.open("time vs no_sim.csv");
     
    for (unsigned int no_sims = no_sims_lower;no_sims<no_sims_upper;no_sims+=no_sims_step){
        myfile<<no_sims<<", ";
        // calculate prices
        myfile<<opt.calculate_price("analytic",s0,r,v,no_sims,false)[0]<<", ";
        myfile<<opt.calculate_price("mc analytic",s0,r,v,no_sims,false)[0]<<", ";
        myfile<<opt.calculate_price("euler",s0,r,v,no_sims,false)[0]<<", ";

        // calculate delta
        myfile<<opt.calculate_delta("analytic",s0,r,v,false)[0]<<", ";
        myfile<<opt.calculate_delta("pw",s0,r,v,no_sims,false)[0]<<", ";
        myfile<<opt.calculate_delta("fd",s0,r,v,no_sims,h,false)[0]<<", ";
        myfile<<opt.calculate_delta("lr",s0,r,v,no_sims,false)[0]<<", ";

        // calculate gamma
        myfile<<opt.calculate_gamma("analytic",s0,r,v,false)[0]<<", ";
        myfile<<opt.calculate_gamma("fd",s0,r,v,no_sims,h,false)[0]<<", ";
        myfile<<opt.calculate_gamma("lr",s0,r,v,no_sims,false)[0]<<", ";
        
        // calculate vega
        myfile<<opt.calculate_vega("analytic",s0,r,v,false)[0]<<", ";
        myfile<<opt.calculate_vega("pw",s0,r,v,no_sims,false)[0]<<", ";
        myfile<<opt.calculate_vega("fd",s0,r,v,no_sims,h,false)[0]<<", ";
        myfile<<opt.calculate_vega("lr",s0,r,v,no_sims,false)[0]<<", ";
        
        // time
        myfile1<<no_sims<<", ";
        // calculate prices
        myfile1<<opt.calculate_price("analytic",s0,r,v,no_sims,false)[1]<<", ";
        myfile1<<opt.calculate_price("mc analytic",s0,r,v,no_sims,false)[1]<<", ";
        myfile1<<opt.calculate_price("euler",s0,r,v,no_sims,false)[1]<<", ";
        
        // calculate delta
        myfile1<<opt.calculate_delta("analytic",s0,r,v,false)[1]<<", ";
        myfile1<<opt.calculate_delta("pw",s0,r,v,no_sims,false)[1]<<", ";
        myfile1<<opt.calculate_delta("fd",s0,r,v,no_sims,h,false)[1]<<", ";
        myfile1<<opt.calculate_delta("lr",s0,r,v,no_sims,false)[1]<<", ";
        
        // calculate gamma
        myfile1<<opt.calculate_gamma("analytic",s0,r,v,false)[1]<<", ";
        myfile1<<opt.calculate_gamma("fd",s0,r,v,no_sims,h,false)[1]<<", ";
        myfile1<<opt.calculate_gamma("lr",s0,r,v,no_sims,false)[1]<<", ";
        
        // calculate vega
        myfile1<<opt.calculate_vega("analytic",s0,r,v,false)[1]<<", ";
        myfile1<<opt.calculate_vega("pw",s0,r,v,no_sims,false)[1]<<", ";
        myfile1<<opt.calculate_vega("fd",s0,r,v,no_sims,h,false)[1]<<", ";
        myfile1<<opt.calculate_vega("lr",s0,r,v,no_sims,false)[1]<<", ";
        
        myfile1<<endl;
        myfile<<endl;
    }
    myfile.close();
    myfile1.close();
    


	return 0;
}
