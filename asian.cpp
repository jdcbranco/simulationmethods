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
#include <fstream>

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

/*************************************************/


class asian_option_geometric{
    double T; // terminal time
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
    
    void MC_euler_pricing(double S0, double r,double v,unsigned int no_sims,vector<double>& res){
        // prameter initlisation
        clock_t c;
        double duration;
        double price = 0;
        double dt = (double)T/N;
        double dt_sqr = pow((double)T/N,0.5); // ADDED: Casting T to double
        double mean_sqr = 0;

        // run the simulation and use antithetic variance reduction
        for(unsigned int i = 0;i < no_sims;i++) {
            vector<double> z = normal_generator(N); //generate normal vector of size N
            double sum_u = 0;
            double sum_d = 0;
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
            mean_sqr += pow(pay_off(exp(sum_u / N)), 2);
            mean_sqr += pow(pay_off(exp(sum_d/ N)), 2);

        }
        // record time
        duration = (clock()-c)/(double)CLOCKS_PER_SEC;
        // return results
        // handle edge cases
        double C = (price / no_sims)*exp(-r * T);
        res.push_back(C); // price
        res.push_back(duration); //time
        res.push_back(C); //mean
        res.push_back((mean_sqr *exp(-r * T) *exp(-r * T) / no_sims) - pow(C, 2)); //variance

    }
    
    void MC_euler_pricing_non_anti(double S0, double r, double v, unsigned int no_sims, vector<double>& res) {
        // prameter initlisation
        clock_t c;
        c = clock();
        double duration;
        double price = 0;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        vector<double> true_val;
        analytic_solution_pricing( S0,  r,  v, true_val);
        double this_price, this_err, err_sum = 0, err_sum_sqr = 0;
        double discount = exp(-r * T);
        
        // run the simulation, NO variance reduction
        for (unsigned int i = 0; i < no_sims; i++) {
            vector<double> z = normal_generator(N); //generate normal vector of size N
            double log_sum = 0;
            double s = S0;
            
            for (unsigned int j = 0;j < N;j++) {
                s += r * s*dt + v * s*dt_sqr*z[j];
                log_sum += log(s);
            }
            this_price = pay_off(exp(log_sum / N)) * discount;
            price += this_price;
            mean_sqr += pow(this_price,2);
            
            this_err = this_price - true_val[0];
            err_sum += this_err;
            err_sum_sqr += pow(this_err,2);
            
        }
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        // return results
        // handle edge cases
        double C = (price / no_sims);
        res.push_back(C); // price
        res.push_back(duration); //time
        res.push_back(err_sum / no_sims ); //mean
        res.push_back((err_sum_sqr / no_sims) - pow(err_sum / no_sims, 2)); //variance
        
    }

    void MC_milstein_pricing(double S0, double r,double v,unsigned int no_sims,vector<double>& res){
        
        // prameter initlisation
        clock_t c;
        c = clock();// start the timer
        double duration;
        double price = 0;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        vector<double> true_val;
        analytic_solution_pricing( S0,  r,  v, true_val);
        double this_price, this_err , err_sum = 0, err_sum_sqr = 0;
        double discount = exp(-r * T);
        
        // run the simulation, NO variance reduction
        for (unsigned int i = 0; i < no_sims; i++) {
            vector<double> z = normal_generator(N); //generate normal vector of size N
            double log_sum = 0;
            double s = S0;
            
            for (unsigned int j = 0;j < N;j++) {
                s += r * s*dt + v * s*dt_sqr*z[j] + 0.5 * v * ( v * s ) * dt * ( z[j]*z[j] - 1 ) ;
                log_sum += log(s);
            }
            this_price = pay_off(exp(log_sum / N)) * discount;
            price += this_price;
            mean_sqr += pow(this_price,2);
            
            this_err = this_price - true_val[0];
            err_sum += this_err;
            err_sum_sqr += pow(this_err,2);
            
        }
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        // return results
        // handle edge cases
        double C = (price / no_sims);
        res.push_back(C); // price
        res.push_back(duration); //time
        res.push_back(err_sum / no_sims ); //mean
        res.push_back((err_sum_sqr / no_sims) - pow(err_sum / no_sims, 2)); //variance
        
    }
    
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
	}
    
    void MC_fd_delta(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {
        
        clock_t c;
        double duration;
        c = clock();// start the timer
        
        double price1_sum = 0, price2_sum = 0, price1, price2, greek_sum = 0, greek_sum_sqr = 0;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        
        double delta, expectation, var;
        vector<double> true_val;
        analytic_solution_delta( S0,  r, v, true_val);
        double this_price_1, this_price_2, this_greek, this_err, err_sum = 0, err_sum_sqr = 0;
        double discount = exp(-r * T);
        
        // run the simulation
        if(icompare(method, "milstein") ) {
            // use milstein if requested
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0;
                double s1 = S0 + h, s2 = S0 - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v * s1*dt_sqr*z[j] + 0.5 * v * ( v * s1 ) * dt * ( z[j]*z[j] - 1 ) ;
                    s2 += r * s2*dt + v * s2*dt_sqr*z[j] + 0.5 * v * ( v * s2 ) * dt * ( z[j]*z[j] - 1 ) ;
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                }
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_greek = (this_price_1 - this_price_2 ) / (2*h) ;

                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
            }
            delta = ( greek_sum / no_sims );
            expectation = err_sum / no_sims ;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else if(icompare(method, "euler")){
            
            // use milstein if requested
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0;
                double s1 = S0 + h, s2 = S0 - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v * s1*dt_sqr*z[j];
                    s2 += r * s2*dt + v * s2*dt_sqr*z[j];
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                }
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_greek = (this_price_1 - this_price_2 ) / (2*h) ;
                
                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
            }
            delta = ( greek_sum / no_sims );
            expectation = err_sum / no_sims ;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else {
            // return ERROR
            cout << "ERROR: scheme method not recognized" << endl;
        }
        
        res.push_back( delta );
        res.push_back( duration );
        res.push_back( expectation );
        res.push_back( var );
    }
    
    void MC_pw_delta(string method, double S0, double r, double v, int no_sims, vector<double> &delta) {
        clock_t c;
        double duration;
        c = clock();// start the timer
        
        // prameter initlisation
        clock_t t;
        double price = 0;
        double dt = (double)T/N;
        double dt_sqr = pow((double)T/N,0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        double sum = 0;
        
        // run the simulation, NO variance reduction
        for (unsigned int i = 0; i < no_sims; i++) {
            vector<double> z = normal_generator(N); //generate normal vector of size N
            double log_sum = log(S0);
            double s = S0;
            
            for (unsigned int j = 0;j < N;j++) {
                s += r * s*dt + v * s*dt_sqr*z[j];
                log_sum += log(s);
            }
            double average = exp(log_sum / (N + 1));
            bool indicator = average > K;
            if( indicator ) {
                sum += average / S0;
                mean_sqr += pow( average / S0,2);
            }
        }
        // record time
        duration = (clock()-t)/(double)CLOCKS_PER_SEC;
        // return results
        // handle edge cases
        delta.push_back((sum/no_sims)); // price
        delta.push_back(duration); //time
        delta.push_back(sum/no_sims); //mean
        delta.push_back((mean_sqr/no_sims)-pow((sum/no_sims),2)); //variance
        
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        
        delta.push_back(0);
        delta.push_back(duration);
        
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
    
    void MC_fd_vega(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {
        
        clock_t c;
        double duration;
        c = clock();// start the timer
        
        double price1_sum = 0, price2_sum = 0, price1, price2, greek_sum = 0, greek_sum_sqr = 0, err_sum = 0, err_sum_sqr = 0;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        
        double vega, expectation, var;
        
        double this_price_1,this_price_2;
        double discount = exp(-r * T);
        double this_err, this_greek;
        vector<double> true_val;
        analytic_solution_vega( S0,  r,  v, true_val);
        
        
        // run the simulation
        if(icompare(method, "milstein") ) {
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0;
                double s1 = S0, s2 = S0;
                double v1 = v + h, v2 = v - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v1 * s1*dt_sqr*z[j] + 0.5 * v1 * ( v1 * s1 ) * dt * ( z[j]*z[j] - 1 ) ;
                    s2 += r * s2*dt + v2 * s2*dt_sqr*z[j] + 0.5 * v2 * ( v2 * s2 ) * dt * ( z[j]*z[j] - 1 ) ;
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                }
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_greek = (this_price_1 - this_price_2 ) / (2*h) ;
                
                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
            }
            vega = ( greek_sum / no_sims );
            expectation = err_sum / no_sims;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else if(icompare(method, "euler")){
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0;
                double s1 = S0, s2 = S0;
                double v1 = v + h, v2 = v - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v1 * s1*dt_sqr*z[j];
                    s2 += r * s2*dt + v2 * s2*dt_sqr*z[j];
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                }
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_greek = (this_price_1 - this_price_2 ) / (2*h) ;
                
                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
            }
            vega = ( greek_sum / no_sims );
            expectation = err_sum / no_sims;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else {
            // return ERROR
            cout << "ERROR: scheme method not recognized" << endl;
        }
        
        res.push_back( vega );
        res.push_back( duration );
        res.push_back( expectation );
        res.push_back( var );
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
    
    void MC_fd_gamma(string method, double S0, double r, double v, int no_sims, vector<double> &res, double h) {
        
        clock_t c;
        double duration;
        c = clock();// start the timer
        
        double greek_sum = 0, greek_sum_sqr = 0, numerator, err_sum = 0, err_sum_sqr;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        
        double res_val, expectation, var;
        
        double this_price_1, this_price_2, this_price_3;
        double discount = exp(-r * T);
        double this_err, this_greek;
        vector<double> true_val;
        analytic_solution_gamma( S0,  r,  v, true_val);
        
        // run the simulation
        if(icompare(method, "milstein") ) {
            // use milstein if requested
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0, log_sum3 = 0;
                double s1 = S0 + h, s2 = S0, s3 = S0 - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v * s1*dt_sqr*z[j] + 0.5 * v * ( v * s1 ) * dt * ( z[j]*z[j] - 1 ) ;
                    s2 += r * s2*dt + v * s2*dt_sqr*z[j] + 0.5 * v * ( v * s2 ) * dt * ( z[j]*z[j] - 1 ) ;
                    s3 += r * s3*dt + v * s3*dt_sqr*z[j] + 0.5 * v * ( v * s3 ) * dt * ( z[j]*z[j] - 1 ) ;
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                    log_sum3 += log(s3);
                }
                
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_price_3 = pay_off(exp(log_sum3 / N))* discount;
                
                numerator = ( this_price_1 - 2*this_price_2 + this_price_3 );
                this_greek = numerator / pow(h,2);
                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
                
            }
            res_val = greek_sum / no_sims;
            expectation = err_sum / no_sims;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else if(icompare(method, "euler")){
            for (unsigned int i = 0; i < no_sims; i++) {
                vector<double> z = normal_generator(N); //generate normal vector of size N, for both price1 AND price2
                double log_sum1 = 0, log_sum2 = 0, log_sum3 = 0;
                double s1 = S0 + h, s2 = S0, s3 = S0 - h;
                
                for (unsigned int j = 0;j < N;j++) {
                    s1 += r * s1*dt + v * s1*dt_sqr*z[j];
                    s2 += r * s2*dt + v * s2*dt_sqr*z[j];
                    s3 += r * s3*dt + v * s3*dt_sqr*z[j];
                    log_sum1 += log(s1);
                    log_sum2 += log(s2);
                    log_sum3 += log(s3);
                }
                
                this_price_1 = pay_off(exp(log_sum1 / N))* discount;
                this_price_2 = pay_off(exp(log_sum2 / N))* discount;
                this_price_3 = pay_off(exp(log_sum3 / N))* discount;
                
                numerator = ( this_price_1 - 2*this_price_2 + this_price_3 );
                this_greek = numerator / pow(h,2);
                greek_sum += this_greek;
                greek_sum_sqr += pow(this_greek,2);
                
                this_err = this_greek - true_val[0];
                err_sum += this_err;
                err_sum_sqr += pow(this_err,2);
                
            }
            res_val = greek_sum / no_sims;
            expectation = err_sum / no_sims;
            var = err_sum_sqr / no_sims - pow( expectation, 2);
            
            // record time
            duration = (clock() - c) / (double)CLOCKS_PER_SEC;
            
        } else {
            // return ERROR
            cout << "ERROR: scheme method not recognized" << endl;
        }
        
        res.push_back( res_val );
        res.push_back( duration );
        res.push_back( expectation );
        res.push_back( var );
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
    vector<double> TEST(double S0, double r,double v,unsigned int no_sims,vector<double>& res){
        
        // prameter initlisation
        clock_t c;
        c = clock();// start the timer
        double duration;
        double price = 0;
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        double mean_sqr = 0;
        
        vector<double> price_vec;
        
        // run the simulation, NO variance reduction
        for (unsigned int i = 0; i < no_sims; i++) {
            vector<double> z = normal_generator(N); //generate normal vector of size N
            double log_sum = 0;
            double s = S0;
            
            for (unsigned int j = 0;j < N;j++) {
                s += r * s*dt + v * s*dt_sqr*z[j] + 0.5 * v * ( v * s ) * dt * ( z[j]*z[j] - 1 ) ;
                log_sum += log(s);
            }
            price_vec.push_back( pay_off(exp(log_sum / N))*exp(-r * T) );
            price += pay_off(exp(log_sum / N));
            mean_sqr += pow(pay_off(exp(log_sum / N )),2);
            
        }
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        // return results
        // handle edge cases
        double C = (price / no_sims)*exp(-r * T);
        res.push_back(C); // price
        res.push_back(duration); //time
        res.push_back(C); //mean
        res.push_back((mean_sqr *exp(-r * T) *exp(-r * T) / no_sims) - pow(C, 2)); //variance
        
        output_file("dist.csv",price_vec);
        return price_vec;
        
    }
    
    
    
    // an asian option has time to maturity T,K,initial price and interest rate
    asian_option_geometric(const string& type,double T,int N, double K):T(T),N(N),K(K){
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
    vector<double> calculate_price(const string& method, double S0, double r, double v, int no_sims = 100000, bool display_flag = true) {
        // return containeer
        vector<double> res;

        //Use selected method to calcuate c_0
        if (icompare(method, "euler")) {
            MC_euler_pricing_non_anti(S0, r, v, no_sims, res);
        } else if (icompare(method,"milstein")){
            MC_milstein_pricing(S0,r,v,no_sims, res);
        } else if (icompare(method,"analytic")){
            analytic_solution_pricing(S0, r, v, res);
        } else {
            cout << "method " << method << " not recognized, please choose another" << endl;
        }
        
        string name = "PRICE";
        if(display_flag) {
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
    vector<double> calculate_delta(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_flag = true ) {
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
        if(display_flag) {
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
    vector<double> calculate_vega(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_flag = true) {
        vector<double> res;
        
        if(icompare(method2, "fd")) {
            MC_fd_vega(method, S0, r, v, no_sims, res, h );
        } else if(icompare(method2, "analytic")) {
            analytic_solution_vega(S0,r,v,res);
        } else {
            // edge case:
        }
        
        string name = "VEGA";
        if(display_flag) {
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
    vector<double> calculate_gamma(string method, string method2, double S0, double r, double v, int no_sims = 100000, double h = 0.01, bool display_flag = true) {
        vector<double> res;
        
        if(icompare(method2, "fd")) {
            MC_fd_gamma(method, S0, r, v, no_sims, res, h );
        } else if(icompare(method2, "analytic")) {
            analytic_solution_gamma(S0,r,v,res);
        } else {
            // edge case:
        }
        
        string name = "GAMMA";
        if(display_flag) {
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

double sum(vector<double> &vec) {
    double sum = 0;
    for(unsigned int i = 0; i < vec.size(); i++ ) {
        sum += vec[i];
    }
    return sum;
}

void normal_convergence_test() {
    int start = 1000, end = 100000, increment = 10000/2;
    vector<double> output;
    for(int m = start; m <= end; m += increment) {
        vector<double> mc_sample;
        for( int n = 0; n < 100; n++ ) {
            vector<double> test_sample = normal_generator(m);
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
    
    /*
     
    opt.calculate_price("analytic",s0,r,v,no_sims);
    opt.calculate_price("euler",s0,r,v,no_sims);
    opt.calculate_price("milstein",s0,r,v,no_sims);
    
    opt.calculate_delta("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_delta("euler","fd",s0,r,v,no_sims, h);
    opt.calculate_delta("milstein","fd",s0,r,v,no_sims, h);
     
    opt.calculate_gamma("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_gamma("euler","fd",s0,r,v,no_sims, h);
    opt.calculate_gamma("milstein","fd",s0,r,v,no_sims, h);
    
    opt.calculate_vega("analytic","analytic",s0,r,v,no_sims, h);
    opt.calculate_vega("euler","fd",s0,r,v,no_sims, h);
    opt.calculate_vega("milstein","fd",s0,r,v,no_sims, h);
    
    */
    
    int start = 10, finish = no_sims, num_of_increments = 10;
    vector<double> no_sims_list = linspace(start,finish,num_of_increments); // gives a series of integers from start to finish, with num_of_increments equal sized steps
    
    price_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list);
    delta_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
    gamma_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
    vega_export_for_changing_M(opt,"milstein",s0,r,v,no_sims_list,h);
    
    cout << "DONE" << endl;
	int dummy;
	cin >> dummy;
	return 0;
}
