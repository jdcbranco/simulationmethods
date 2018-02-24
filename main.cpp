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
    double T;
    double sigma;

public:
    Derivatives(double strike = 100.0, double initial_price = 100.0,double T =1.0 , double sigma = 0.4) {
        K = strike;
        S0 = initial_price;
        this -> T = T;
        this -> sigma = sigma ;
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

public :
    model (Derivatives d,double r) {
        this ->d = d;
        this->r = r;
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

    double mean_vector(vector<double> vect) {
        double sum = 0.0 ;
        for (unsigned int i = 0 ; i<vect.size();i=i+1) {
            sum = sum+vect[i] ;
        }
        return sum/vect.size();
    }

    double variance_vector(vector<double> v){
    		int n = v.size();
    		double mean=mean_vector(v);
    		double sum=0.0;
    		for(int i=0;i<n;i++){
    			sum+=pow(v[i]-mean,2);
    		}
    		sum/=(n-1);
    		return pow(sum,0.5);
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

public :

    MonteCarlo (Derivatives d,double r, int N) : model(d,r) {
        this -> N = N;
    }

    // generate the price path

   /*
       void euler_path() {
        int i;
        price_path = {d.get_S0()} ;
        vector<double> normal_random = generate_normal(0.0, 1.0, N) ;
        for (i=1 ; i <N ; i = i+1) {
            double S_i = price_path[i-1] + r*price_path[i-1]*d.get_T()/N + d.get_sigma()*price_path[i-1]*pow(d.get_T()/N , 0.5)*normal_random[i];
            price_path.push_back(S_i);
        }
    }
    */


    vector<double> payoff_call(Derivatives d, int n_simulations,int epsilon=1) {  // Epsilon will be equal to +1 or -1 : this is to use the antithetic method for variance
        vector<double> payoff ;
        int i ;
        vector<double> normal_random = generate_normal(0.0, 1.0, n_simulations) ;
        for (i = 0 ; i < n_simulations ; i = i+1) {
            double S_T = d.get_S0()*exp( (r- pow(d.get_sigma(),2.0)/2)*d.get_T() + epsilon*d.get_sigma()*pow(d.get_T(), 0.5)*normal_random[i]);
            if (S_T-d.get_K() <=0 ) {
                payoff.push_back(0.0)  ;
            }
            else {
                payoff.push_back(exp(-r*d.get_T())*(S_T-d.get_K()));
            }
        }
        return payoff;
    }

    double CalculVariance(Derivatives d, int n_simulations){
    		vector<double> payoffs = payoff_call(d,n_simulations);
    		double variance= variance_vector(payoffs);
    		return variance;

    }


    double CalculPrice(Derivatives d, int n_simulations, int epsilon =1) {
        vector<double> payoffs = payoff_call(d,n_simulations, epsilon);
        double mean_payoffs = mean_vector(payoffs);
        return mean_payoffs;
    }


    double CalculDelta(Derivatives d, int n_simulations,double h) {
        // same derivative with price S0-h
        Derivatives d1 ;
        d1 = d;
        d1.set_S0(d.get_S0() -h) ;
        // same derivative with price S0+h
        Derivatives d2;
        d2 = d ;
        d2.set_S0(d.get_S0() +h) ;

        // payoffs of derivatives with initial price S0-h and S0+h
        vector<double> payoff1 = payoff_call(d1, n_simulations) ;
        vector<double> payoff2 = payoff_call(d2, n_simulations) ;

        double mean_payoff1 = mean_vector(payoff1); // mean of payoffs of the derivative with initial price S0-h
        double mean_payoff2 = mean_vector(payoff2); // mean of payoffs of the derivative with initial price S0+h
        return ((mean_payoff2-mean_payoff1)/(2*h));

    }

    double CalculDelta_2(Derivatives d, int n_simulations){

    	vector<double> payoff1 = payoff_call(d, n_simulations);
    	    	double K=d.get_K();
    	    	double S0=d.get_S0();
    	    	double sum=0;
    	    	for (unsigned int i=0;i<payoff1.size();i++){
    	    		if(payoff1[i]>0){
    	    			sum+=(payoff1[i]+K)/S0;
    	    		}
    	    	}
    	    	return (sum*exp(-r*d.get_T()))/n_simulations;


    }

    double CalculGamma(Derivatives d, int n_simulations,double h) {
        // same derivative with price S0-h
        Derivatives d1 ;
        d1 = d;
        d1.set_S0(d.get_S0()-h) ;
        // same derivative with price S0+h
        Derivatives d2;
        d2 = d ;
        d2.set_S0(d.get_S0()+h) ;

        // payoffs of derivatives with initial price S0-h, S0 and S0+h
        vector<double> payoff1 = payoff_call( d1, n_simulations);
        vector<double> payoff2 = payoff_call( d, n_simulations) ;
        vector<double> payoff3 = payoff_call( d2 , n_simulations) ;

        // mean of payoffs
        double mean_payoff1 = mean_vector(payoff1); // initial price S0-h
        double mean_payoff2 = mean_vector(payoff2); // initial price S0
        double mean_payoff3 = mean_vector(payoff3); // initial price S0+h

        return ((-mean_payoff1+mean_payoff3-2*mean_payoff2)/(h*h));

    }


    double CalculVega(Derivatives d,int n_simulations,double h) {
        // same derivative with volatility sigma-h
        Derivatives d1 ;
        d1 = d;
        d1.set_sigma(d.get_sigma() -h) ;
        // same derivative with volatility sigma+h
        Derivatives d2;
        d2 = d ;
        d2.set_sigma(d.get_sigma() +h) ;

        // payoffs of derivatives with initial volatility sigma-h and sigma+h
        vector<double> payoff1 = payoff_call(d1, n_simulations) ;
        vector<double> payoff2 = payoff_call(d2, n_simulations) ;

        double mean_payoff1 = mean_vector(payoff1); // mean of payoffs of the derivative with initial volatility sigma-h
        double mean_payoff2 = mean_vector(payoff2); // mean of payoffs of the derivative with initial volatility sigma+h
        return ((mean_payoff2-mean_payoff1)/(2*h));
    }


    double CalculVega_2(Derivatives d,int n_simulations){
    	vector<double> payoff1 = payoff_call(d, n_simulations);
    	        	double K=d.get_K();
    	        	double S0=d.get_S0();
    	        	double sum=0;
    	        	double sigma=d.get_sigma();
    	        	double ST;
    	        	double T=d.get_T();
    	        	for (unsigned int i=0;i<payoff1.size();i++){
    	        		if(payoff1[i]>0){
    	        			ST=payoff1[i]+K;
    	        			sum+=((ST)/sigma)*(log(ST/S0)-(r+0.5*pow(sigma,2)*T));
    	        		}
    	        	}
   	return (sum*exp(-r*d.get_T()))/n_simulations;
    }


};


class Black_Scholes :public model{

public:
    Black_Scholes(Derivatives d,double r) : model(d,r) { }

    double normalCDF(double value) {
        return 0.5 * erfc(-value / sqrt(2));
    }

    double CalculPrice() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
        double d2 = d1 - d.get_sigma() * sqrt(d.get_T());
        return (d.get_S0()*normalCDF(d1) - normalCDF(d2)*d.get_K()*exp(-r*d.get_T()) );
    }

    double CalculDelta() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
        return (normalCDF(d1));
    }

    double CalculGamma() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
        return (exp(-d1*d1/2)/sqrt(2*M_PI) /(d.get_S0()*d.get_sigma()*sqrt(d.get_T())) ) ;
    }

    double CalculVega() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + d.get_sigma()*d.get_sigma() / 2)*d.get_T()) / (d.get_sigma()*sqrt(d.get_T()));
        return(exp(-d1*d1/2)/sqrt(2*M_PI) * d.get_S0() *sqrt(d.get_T())) ;
    }
};

class MC_Asian_Call : public model {
	int N; // number of path incrememnts, to be input as a model parameter

	// TEMP FIX:
	double sigma = d.get_sigma();
	double T = d.get_T();
public :
	MC_Asian_Call(Derivatives d, double r, int N) : model(d, r) {
		// constructor to include N in additional to existing 
		this -> N = N;

		// TEMP FIX for the merge
		sigma = d.get_sigma();
		T = d.get_T();
		
	}

	double geo_average(vector<double> &vec) {
		// function to caluclate geometric average, given vector, uses finds the average of the log values first, then gives exp( average(log values) )
		double log_sum = 0.0; 
		for (unsigned int i = 0; i < vec.size(); i++) {
			log_sum = log_sum + log(vec[i]);
		}
		return exp(log_sum / vec.size());
	}

	vector<double> generate_path(Derivatives d) {
		double h = T / N;
		vector<double> path(N+1);
		vector<double> Z = generate_normal(0.0, 1.0, N);

		path[0] = d.get_S0();

		for (int i = 0; i < N; i++) {
			path[i+1] = path[i] + r * path[i]*h + sigma * path[i]*sqrt(h)*Z[i];
		}
		return path;
	}

	double CalculPrice(int M, Derivatives d) {
		// M is number of MC simulations;

		// payoff for asian option is: A_T = ( geo_average( price_path ) - K )^+

		vector<double> A_T_vec;
		for (int i = 0; i < M; i++) {
			// M times of MC simulations
			vector<double> price_path = generate_path(d);
			double average = geo_average(price_path);		
			A_T_vec.push_back(((average - d.get_K()) >= 0) * (average - d.get_K()));
		}

		return mean_vector(A_T_vec);
	}

	double CalculPrice(int M) {
		return CalculPrice(M, d);
	}

	double CalculDelta(int M, double h) {
		// same derivative with stock price S0-h
		Derivatives d1;
		d1 = d;
		d1.set_S0(d.get_S0() - h);
		// same derivative with stock price S0+h
		Derivatives d2;
		d2 = d;
		d2.set_S0(d.get_S0() + h);

		// option prices with initial stock price S0-h and S0+h
		double price1 = CalculPrice(M, d1);
		double price2 = CalculPrice(M, d2);
		return ((price2 - price1) / (2 * h));

	}

	double CalculDelta_2(int M) {
		double K = d.get_K();
		double S0 = d.get_S0();
		double sum = 0.0;
		for (int i = 0;i<M;i++) {
			vector<double> path = generate_path(d); 
			double average = geo_average(path);
			sum += (average > K) * average / S0;
		}
		return exp(-r * T)* sum / M;
	}

	double CalculDelta_3(int M) {
		double K = d.get_K();
		double S0 = d.get_S0();
		double sum = 0.0;
		double t = 0.0;
		double mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
		double sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);

		for (int i = 0; i < M; i++) {
			vector<double> path = generate_path(d);
			double average = geo_average(path);
			sum += (average > K) * (average - K ) * ( log( average / S0 ) - mu_bar );
		}

		return exp(-r * T) * sum / M / ( S0* pow(sigma_bar,2) );
	}

};


class CF_Asian_Call : public model {
	double t = 0.0;
	double G_t, mu_bar, sigma_bar, d1,d2;

	// TEMP FIX:
	double sigma, T;

public:
	CF_Asian_Call(Derivatives d, double r) : model(d, r) {
		
		// TEMP FIX:
		sigma = d.get_sigma();
		T = d.get_T();
		
		
		G_t = d.get_S0();
		mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
		sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);
		d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(d.get_S0()) + mu_bar - log(d.get_K()));
		d1 = d2 + sigma_bar;
	}
	double normalCDF(double value) {
		// PROPOSAL: move this higher up as a public function instead of a member function?
		return 0.5 * erfc(-value / sqrt(2));
	}

	double normalPDF(double value) {
		return (1 / sqrt(2 * M_PI)) * exp(-0.5 * pow(value, 2));
	}

	double CalculPrice() {
		return exp(-r * (T - t)) * (pow(G_t, t / T) * pow(d.get_S0(), (T - t) / T) * exp(mu_bar + pow(sigma_bar, 2) / 2) * normalCDF( d1 ) - d.get_K()*normalCDF(d2));
	}

	double CalculDelta() {
		return exp(mu_bar + pow(sigma_bar, 2) / 2)* normalCDF(d1) + exp(mu_bar + pow(sigma_bar, 2) / 2) / sigma_bar * normalPDF(d1) - d.get_K() / (sigma_bar * d.get_S0()) * normalPDF(d2); 
	}
};

int main() {
	double strike=100;
	int number_simulations=1000000;
	double sigma=0.4;
	double r=0.05;
    Derivatives call(strike,100,1.0,sigma);
    MonteCarlo mc(call ,r,100) ;
    Black_Scholes bs(call ,r) ;
    // price of the option
    double price_call_mc = mc.CalculPrice(call, number_simulations);
    cout << "The Monte_Carlo price of the call is equal to "<< price_call_mc << endl;
    double price_call_bs = bs.CalculPrice();
    cout << "The BS price of the call is equal to "<< price_call_bs << endl;

    double price_call_MC_antit=(mc.CalculPrice(call, number_simulations)+mc.CalculPrice(call, number_simulations,-1))/2;
    cout << "The MC price of the call using antithetic variables is equal to "<< price_call_MC_antit << endl;

    vector<double> v=mc.payoff_call(call, number_simulations,+1);
    vector<double> vminus=mc.payoff_call(call, number_simulations,-1);
    for(int i=0;i<50;i++){
    		cout<<v[i]<<" ";

    }
    cout<<endl;

    for(int i=0;i<50;i++){
        		cout<<vminus[i]<<" ";

        }
    vector<double> v_minus=mc.payoff_call(call, number_simulations,-1);
    v.insert(v.end(),v_minus.begin(),v_minus.end());

    double variance_antiT=mc.variance_vector(v);
    cout << "The Monte_Carlo variance of the call using Antitetic variables is equal to "<< variance_antiT << endl;


    //mean
    double variance_MC=mc.CalculVariance(call,number_simulations);
    cout << "The Monte_Carlo variance of the call is equal to "<< variance_MC << endl;


    // delta
    /*double delta_mc = mc.CalculDelta(call, 1000000, 0.1);
    cout << "The Monte_Carlo delta is equal to "<< delta_mc << endl;*/
    double delta_bs = bs.CalculDelta();
    cout << "The BS delta is equal to "<< delta_bs << endl;
    double delta_new=mc.CalculDelta_2(call,number_simulations);
    cout << "The Monte_Carlo delta_1 is equal to "<< delta_new << endl;

    // gamma
    double gamma_mc = mc.CalculGamma(call, number_simulations, 0.001);
    cout << "The MonteCarlo gamma is equal to "<< gamma_mc << endl;
    double gamma_bs = bs.CalculGamma();
    cout << "The BS gamma is equal to "<< gamma_bs << endl;

    // vega
    double vega_mc = mc.CalculVega_2(call, number_simulations);
    cout << "The Monte Carlo vega is equal to " << vega_mc << endl;
    double vega_bs = bs.CalculVega();
	cout << "The BS vega is equal to " << vega_bs << endl;


	// Asian options
	cout << endl;
	for (int i = 0; i < 20;i++) {
		cout << "-";
	}
	cout << endl;
	cout << "Asian Options" << endl;
	for (int i = 0; i < 20;i++) {
		cout << "-";
	}

	Derivatives asian_call(100.0, 100.0, 1.0, 0.4);
	MC_Asian_Call mc_asian_call(asian_call, 0.05, 100);
	CF_Asian_Call cf_asian_call(asian_call, 0.05);

	// price of the option
	cout << endl;
	double asian_call_price_mc = mc_asian_call.CalculPrice(100000);
	double asian_call_price_cf = cf_asian_call.CalculPrice();
	cout << "The Monte Carlo price of the asian call is equal to " << asian_call_price_mc << endl;
	cout << "The Closed form price of the asian call is equal to " << asian_call_price_cf << endl;

	cout << endl;
	double asian_call_delta_mc;
	// double asian_call_delta_mc = mc_asian_call.CalculDelta(10000, 0.1);
	asian_call_delta_mc = mc_asian_call.CalculDelta_2(100000);
	cout << "The Monte Carlo delta of the asian call is equal to " << asian_call_delta_mc << endl;
	asian_call_delta_mc = mc_asian_call.CalculDelta_3(100000);
	cout << "The Monte Carlo delta of the asian call is equal to " << asian_call_delta_mc << endl;
	double asian_call_delta_cf = cf_asian_call.CalculDelta();	
	cout << "The Closed form delta of the asian call is equal to " << asian_call_delta_cf << endl;


	int dummy;
	cin >> dummy;
    return 0;
}
 