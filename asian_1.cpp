#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "Normal.h"

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

/**************************************************************************************************/

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

    void analytic_solution_pricing(double S0, double r, double v, vector<double>& res){
        double mu_a = (double)T*(r-v*v/2.0)*((N+1.0)/2.0*N);
        double var_a = (double)T*v*v*(1.0/N+(N-1.0)*(2.0*N-1)/(6.0*N*N));
        double sd_a = sqrt(var_a);
        double d2 = (1.0/sd_a)*(log(S0/K)+mu_a);
        double d1 = d2+sd_a;
        double price = S0 * exp(mu_a + 0.5*var_a)*normalCDF(d1) - K*normalCDF(d2);

        // return results
        res.push_back(price);
        res.push_back(0);
    }

    void MC_fd_delta(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {
        clock_t c;
        double duration;
        double mu_a = T*(r-v*v/2)*((N+1)/2*N);
        double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
        double sd_a = sqrt(var_a);
        double discount = exp(-r * T);
        double price1,price2 = 0;
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
        price1 = discount*price2/no_sims;
        delta = (price1+price2)/(2*h);

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;

        res.push_back( delta );
        res.push_back( duration );

    }

    void MC_pw_delta(string method, double S0, double r, double v, int no_sims, vector<double> &res) {
        clock_t c;
        double duration;
        double mu_a = T*(r-v*v/2)*((N+1)/2*N);
        double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
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

    void MC_lr_delta(double S0, double r, double v, int no_sims, vector<double> &res){
        clock_t c;
        double duration;
        double mu_a = T*(r-v*v/2)*((N+1)/2*N);
        double sd_a = sqrt(T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N)));
        double buffer = 0;
        double mean = 0;
        double mean_sqr = 0;

        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(int i=0;i<no_sims;i++){
          buffer = pay_off(S0 * exp(mu_a + sd_a*z[i]))*z[i]/(sd_a*S0);
          mean += buffer;
          mean_sqr += buffer * buffer;
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        mean /= no_sims;
        mean_sqr /= no_sims;

        // return the display_results
        res.push_back(mean*exp(-r*T)); //delta
        res.push_back(duration); //time
        res.push_back(mean); //mean
        res.push_back(mean_sqr - mean * mean); //variance


    }

    void analytic_solution_delta(double S0, double r,double sigma, vector<double> &delta) {
        clock_t c;
        double duration;
        c = clock(); // start the timer

        double t = 0;

        double G_t = S0;
        double mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
        double sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);
        double d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(S0) + mu_bar - log(K));
        double d1 = d2 + sigma_bar;

        double delta_val = exp(mu_bar + pow(sigma_bar, 2) / 2)* normalCDF(d1) + exp(mu_bar + pow(sigma_bar, 2) / 2) / sigma_bar * normalPDF(d1) - K / (sigma_bar * S0) * normalPDF(d2);

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;

        delta.push_back( delta_val );
        delta.push_back( duration );
    }

    void MC_fd_gamma(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {

      clock_t c;
      double duration;
      double mu_a = T*(r-v*v/2)*((N+1)/2*N);
      double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
      double sd_a = sqrt(var_a);
      double discount = exp(-r * T);
      double price1,price2,price3 = 0;
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
      price3 = discount*price2/no_sims;
      gamma = (price1-2*price2+price3)/(h*h);

      // record time
      duration = (clock() - c) / (double)CLOCKS_PER_SEC;

      res.push_back( gamma );
      res.push_back( duration );
    }

    void MC_lr_gamma(double S0, double r, double v, int no_sims, vector<double> &res){
        clock_t c;
        double duration;
        double mu_a = T*(r-v*v/2)*((N+1)/2*N);
        double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
        double sd_a = sqrt(T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N)));
        double S0_sqr = S0*S0;
        double buffer = 0;
        double mean = 0;
        double mean_sqr = 0;


        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(int i=0;i<no_sims;i++){
          buffer = pay_off(S0 * exp(mu_a + sd_a*z[i]))*((z[i]*z[i]-1)/(S0_sqr*var_a)-z[i]/(S0_sqr*sd_a));
          mean += buffer;
          mean_sqr += buffer * buffer;
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        mean /= no_sims;
        mean_sqr /= no_sims;


        // return the display_results
        res.push_back(mean*exp(-r*T)); //gamma
        res.push_back(duration); //time
        res.push_back(mean); //mean
        res.push_back(mean_sqr - mean * mean); //variance
    }

    void analytic_solution_gamma(double S0, double r, double sigma, vector<double> &res) {
        clock_t c;
        double duration;
        c = clock(); // start the timer

        double t = 0;

        double G_t = S0;
        double mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
        double sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);
        double d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(S0) + mu_bar - log(K));
        double d1 = d2 + sigma_bar;

        double d_d = 1.0 / ( sigma * S0 );

        double res_val =    exp( mu_bar + pow( sigma_bar,2 ) ) * normalPDF( d1 ) * d_d
                            + exp( mu_bar + pow(sigma_bar,2)) / sigma_bar  * ( -d1 ) * normalPDF(d1) * d_d
                            - ( -d2 * d_d * K * normalPDF(d2) * sigma_bar * S0 - K * normalPDF(d2) * sigma_bar ) / pow( sigma_bar * S0, 2);

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;

        res.push_back( res_val );
        res.push_back( duration );

    }

    void MC_lr_vega(double S0, double r, double v, int no_sims, vector<double> &res){
        clock_t c;
        double duration;
        double mu_a = T*(r-v*v/2)*((N+1)/2*N);
        double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
        double sd_a = sqrt(T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N)));
        double g =  T*v*(1/N+(N-1)*(2*N-1)/(6*N*N))/sd_a;
        double buffer = 0;
        double mean = 0;
        double mean_sqr = 0;


        // start the timer
        c = clock();

        // run simulation
        vector<double> z = normal.generate(no_sims);
        for(int i=0;i<no_sims;i++){
          buffer = pay_off(S0 * exp(mu_a + sd_a*z[i]))*(g*(pow(z[i],2)-1)/sd_a);
          mean += buffer;
          mean_sqr += buffer * buffer;
        }

        // stop the timer
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;

        // calculate the results
        mean /= no_sims;
        mean_sqr /= no_sims;

        // return the display_results
        res.push_back(mean*exp(-r*T)); //vega
        res.push_back(duration); //time
        res.push_back(mean); //mean
        res.push_back(mean_sqr - mean * mean); //variance
    }

    void MC_pw_vega(double S0, double r, double v, int no_sims, vector<double> &res){
      clock_t c;
      double duration;
      double mu_a = T*(r-v*v/2)*((N+1)/2*N);
      double var_a = T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N));
      double sd_a = sqrt(T*v*v*(1/N+(N-1)*(2*N-1)/(6*N*N)));
      double g =  T*v*(1/N+(N-1)*(2*N-1)/(6*N*N))/sd_a;
      double A = 0;
      double buffer = 0;
      double mean = 0;
      double mean_sqr = 0;


      // start the timer
      c = clock();

      // run simulation
      vector<double> z = normal.generate(no_sims);
      for(int i=0;i<no_sims;i++){
        // calculate the average
        A = S0 * exp(mu_a + sd_a*z[i]);
        // if the average is too small
        if(A<=K){
          continue;
        }

        buffer = z[i]*exp(mu_a + sd_a *z[i])*g;
        mean += buffer;
        mean_sqr += buffer * buffer;
      }

      // stop the timer
      duration = (clock()-c)/(double)CLOCKS_PER_SEC;

      // calculate the results
      mean /= no_sims;
      mean_sqr /= no_sims;

      // return the display_results
      res.push_back(mean*exp(-r*T)); //vega
      res.push_back(duration); //time
      res.push_back(mean); //mean
      res.push_back(mean_sqr - mean * mean); //variance
    }

    void MC_fd_vega(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {
      clock_t c;
      double duration;
      double vega = 0;
      // prameters for the first price
      double mu_a1 = T*(r-(v+h)*(v+h)/2)*((N+1)/2*N);
      double var_a1 = T*(v+h)*(v+h)*(1/N+(N-1)*(2*N-1)/(6*N*N));
      double sd_a1 = sqrt(var_a1);

      // prameters for the second price
      double mu_a2 = T*(r-(v-h)*(v-h)/2)*((N+1)/2*N);
      double var_a2 = T*(v-h)*(v-h)*(1/N+(N-1)*(2*N-1)/(6*N*N));
      double sd_a2 = sqrt(var_a2);
      double discount = exp(-r * T);
      double price1,price2 = 0;

      c = clock();// start the timer

      // run simulation
      vector<double> z = normal.generate(no_sims);
      for (unsigned int i = 0; i < no_sims; i++) {
        price1 += pay_off(S0 * exp(mu_a1 + sd_a1 *z[i]));
        price2 += pay_off(S0 * exp(mu_a2 + sd_a2 *z[i]));
      }

      // calculate detla
      price1 = discount*price1/no_sims;
      price1 = discount*price2/no_sims;
      vega = (price1+price2)/(2*h);

      // record time
      duration = (clock() - c) / (double)CLOCKS_PER_SEC;

      res.push_back( vega );
      res.push_back( duration );
    }

    void analytic_solution_vega(double S0, double r, double sigma, vector<double> &vega) {
        clock_t c;
        double duration;
        c = clock(); // start the timer

        double t = 0;

        double G_t = S0;
        double mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
        double sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);
        double d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(S0) + mu_bar - log(K));
        double d1 = d2 + sigma_bar;

        double d_mu_bar = -sigma*T / 2;
        double d_sigma_bar = sqrt(T/3);

        double d_d2 = ( d_mu_bar * sigma_bar - ( log(S0) + mu_bar - log(K) ) * d_sigma_bar ) / pow( sigma_bar, 2);
        double d_d1 = d_d2 + d_sigma_bar;

        double vega_val =  S0 * exp( mu_bar + pow(sigma_bar,2) / 2 ) * ( d_mu_bar + sigma_bar * d_sigma_bar ) * normalCDF( d1 )
                            + S0 * exp( mu_bar + pow(sigma_bar,2) / 2 ) * normalPDF( d1 ) * d_d1
                            - K * normalPDF( d2 ) * d_d2;

        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;

        vega.push_back( vega_val );
        vega.push_back( duration );

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

public :

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
    vector<double> calculate_price(const string& method, double S0, double r, double v, int no_sims = 100000, bool display_results = true) {
        // return containeer
        vector<double> res;

        //Use selected method to calcuate c_0
        if (icompare(method, "euler")) {
            MC_euler_pricing(S0, r, v, no_sims, res);
        } else if (icompare(method,"analytic")){
            analytic_solution_pricing(S0, r, v, res);
        } else {
            cout << "method " << method << " not recognized, please choose another" << endl;
        }

        string name = "PRICE";
        if(display_results) {
            display_opening(name,method);

            display_section_split();

            cout << "   S0       = " << S0 << endl;
            cout << "   r        = " << r  << endl;
            cout << "   v        = " << v  << endl;

            if( !icompare(method,"analytic") ) {
                cout << "   M        = " << no_sims << endl;
            }
            display_section_split();
            res_print(res);
            display_ending();
        }

        return res;

    }

    // this function calculate delta at time 0 and return analytic statistics
    vector<double> calculate_delta(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_results = true ) {
        // method is between euler/ emilistein
        // method2 is the user choice of between finite difference/ pathwise/ likelihood ratio
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_delta(method, S0, r, v, no_sims, res, h);
            //MC_fd_delta_alt(method, S0, r, v, no_sims, delta, h);


        } else if(icompare(method2, "pw")) {
            MC_pw_delta(method, S0, r, v, no_sims, res);


        } else if(icompare(method2, "analytic") ) {
            analytic_solution_delta(S0,r,v,res);
        } else {
            // edge case:
        }

        string name = "DELTA";
        if(display_results) {
            display_opening(name,method,method2);

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
    vector<double> calculate_vega(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_results = true) {
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_vega(method, S0, r, v, no_sims, res, h );
        } else if(icompare(method2, "analytic")) {
            analytic_solution_vega(S0,r,v,res);
        } else {
            // edge case:
        }

        string name = "VEGA";
        if(display_results) {
            display_opening(name,method,method2);

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
    vector<double> calculate_gamma(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_results = true) {
        vector<double> res;

        if(icompare(method2, "fd")) {
            MC_fd_gamma(method, S0, r, v, no_sims, res, h );
        } else if(icompare(method2, "analytic")) {
            analytic_solution_gamma(S0,r,v,res);
        } else {
            // edge case:
        }

        string name = "GAMMA";
        if(display_results) {
            display_opening(name,method,method2);

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













void price_export_for_changing_M( asian_option_geometric &opt, string method, double s0, double r, double v, vector<double> &no_sims_list ) {
    vector< vector<double> > output;
    vector<double> res;
    int M;
    for (unsigned int i = 0; i<no_sims_list.size(); i++ ) {
        M = no_sims_list[i];

        res = opt.calculate_price(method, s0, r,v , M, false);
        res.push_back( M );
        output.push_back( res );

        cout << "Price stats: " <<  (int) ( (double)i/no_sims_list.size() *100 ) << "% DONE" << endl;
    }
    string name = "price.csv";
    output_file(name,output);
    cout << "Price stats output as: " << name << endl << endl;
}

void delta_export_for_changing_M(asian_option_geometric &opt, string method, double s0, double r, double v, vector<double> &no_sims_list, double h = 0.01) {
    vector< vector<double> > output;
    vector<double> res;
    int M;
    for (unsigned int i = 0; i<no_sims_list.size(); i++ ) {
        M = no_sims_list[i];

        res = opt.calculate_delta(method, "fd", s0, r, v, M, h, false);
        res.push_back( M );
        output.push_back( res );

        cout << "Delta stats: "<<  (int) ( (double)i/no_sims_list.size() *100 ) << "% DONE" << endl;
    }
    string name = "delta.csv";
    output_file(name,output);
    cout << "Delta stats output as: " << name << endl << endl;
}

void gamma_export_for_changing_M(asian_option_geometric &opt, string method, double s0, double r, double v, vector<double> &no_sims_list, double h = 0.01) {
    vector< vector<double> > output;
    vector<double> res;
    int M;
    for (unsigned int i = 0; i<no_sims_list.size(); i++ ) {
        M = no_sims_list[i];

        res = opt.calculate_gamma(method, "fd", s0, r,v , M, h, false);
        res.push_back( M );
        output.push_back( res );

        cout << "gamma stats: "<<  (int) ( (double)i/no_sims_list.size() *100 ) << "% DONE" << endl;
    }
    string name = "gamma.csv";
    output_file(name,output);
    cout << "gamma stats output as: " << name << endl << endl;
}

void vega_export_for_changing_M(asian_option_geometric &opt, string method, double s0, double r, double v, vector<double> &no_sims_list, double h = 0.01) {
    vector< vector<double> > output;
    vector<double> res;
    int M;
    for (unsigned int i = 0; i<no_sims_list.size(); i++ ) {
        M = no_sims_list[i];

        res = opt.calculate_vega(method, "fd", s0, r,v , M, h, false);
        res.push_back( M );
        output.push_back( res );

        cout << "vega stats: "<<  (int) ( (double)i/no_sims_list.size() *100 ) << "% DONE" << endl;
    }
    string name = "vega.csv";
    output_file(name,output);
    cout << "vega stats output as: " << name << endl << endl;
}


double mean(vector<double> &vec ) {
    double sum = 0;
    for(unsigned int i = 0; i < vec.size(); i++ ) {
        sum += vec[i];
    }
    return sum / vec.size();
}

double var(vector<double> &vec) {
    double sum_sqr = 0;
    double m = mean(vec);
    for(unsigned int i = 0; i < vec.size(); i++ ){
        sum_sqr += pow( vec[i] , 2);
    }
    return sum_sqr / vec.size() - pow(m,2);
}

void normal_convergence_test() {
    unsigned int start = 1000, end = 100000, increment = 10000/2;
    vector<double> output;
    for(int m = start; m <= end; m += increment) {
        vector<double> mc_sample;
        for( int n = 0; n < 100; n++ ) {
            vector<double> test_sample = normal.generate(m);
            mc_sample.push_back( mean(test_sample) );
        }

        // cout << "mean = " << mean(test_sample) << endl;
        cout << "var  = " << var(mc_sample) << endl;
        output.push_back( var(mc_sample) );
        // no_sims_list.push_back(m);
    }

    output_file("test.csv",output);

}

vector<double> linspace( int start = 10, int end = 100000, int num_incre = 10) {
    vector<double> output;
    for(int m = start; m <= end; m += (end-start)/num_incre ) {
        output.push_back(m);
    }
    return output;
}


int main() {

//    This code helps identify any problem with the normal random number generator
//    vector<double> lst1 = normal.generate(10);
//    vector<double> lst2 = normal.generate(10);
//    copy(lst1.begin(), lst1.end(), ostream_iterator<double>(cout," "));
//    cout << endl;
//    copy(lst2.begin(), lst2.end(), ostream_iterator<double>(cout," "));
//    cout << endl;
//    return 0;

	int T = 1;
	unsigned int N = 100;
	double K = 100;
	string method = "euler";
	string type = "call";
	double s0 = 100;
	double r = 0.05;
	double v = 0.4;
	int no_sims = 100000;
    double h = 0.01;

    // cout << "Input number of simulation: ";
    // cin >> no_sims;


	srand(time(NULL));
	asian_option_geometric opt(type,T,N,K);

    opt.calculate_price("analytic",s0,r,v,no_sims);
    opt.calculate_price("euler",s0,r,v,no_sims);

    /*
    
    opt.calculate_delta("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_delta("euler","fd",s0,r,v,no_sims, h);

    opt.calculate_gamma("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_gamma("euler","fd",s0,r,v,no_sims, h);

    opt.calculate_vega("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_vega("euler","fd",s0,r,v,no_sims, h);

    
    /*
    int start = 10, finish = no_sims, num_of_increments = 10;
    vector<double> no_sims_list = linspace(start,finish,num_of_increments); // gives a series of integers from start to finish, with num_of_increments equal sized steps

    price_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list);
    delta_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
    gamma_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
    vega_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
     
     */
     
    cout << "DONE" << endl;
	int dummy;
	cin >> dummy;
	return 0;
}
