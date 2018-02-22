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
    Derivatives(double strike = 100.0, double initial_price = 100.0) {
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
    Derivatives d;
    double r;
    double T;
    double sigma;

public :
    model (Derivatives d,double r,double T,double sigma) {
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
    		return pow(sum/double(n-1),0.5);
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

    MonteCarlo (Derivatives d,double r,double T,double sigma, int N) : model(d,r,T,sigma) {
        this -> N = N;
    }

    // generate the price path

   /* void euler_path() {
        int i;
        vector<double> normal_random = generate_normal(0.0, 1.0, N) ;
        for (i=1 ; i <N ; i = i+1) {
            double S_i = price_path[i-1] + r*price_path[i-1]*T/N + sigma*price_path[i-1]*pow(T/N , 0.5)*normal_random[i];
            price_path.push_back(S_i);
        }
    }*/


    vector<double> payoff_call(Derivatives d, int n_simulations) {
        vector<double> payoff ;
        int i ;
        vector<double> normal_random = generate_normal(0.0, 1.0, n_simulations) ;
        for (i = 0 ; i < n_simulations ; i = i+1) {
            double S_T = d.get_S0()*exp( (r- pow(sigma,2.0)/2)*T + sigma*pow(T, 0.5)*normal_random[i]);
            if (S_T-d.get_K() < 0 ) {
                payoff.push_back(0.0)  ;
            }
            else {
                payoff.push_back(exp(-r*T)*(S_T-d.get_K()));
            }
        }
        return payoff;
    }

    double CalculVariance(Derivatives d, int n_simulations){
    		vector<double> payoffs = payoff_call(d,n_simulations);
    		double variance= variance_vector(payoffs);
    		return variance;

    }


    double CalculPrice(Derivatives d, int n_simulations) {
        vector<double> payoffs = payoff_call(d,n_simulations);
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
    	return (sum*exp(-r*T))/n_simulations;



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

        return ((mean_payoff1+mean_payoff3-2*mean_payoff2)/(h*h));

    }





};


class Black_Scholes :public model{

public:
    Black_Scholes(Derivatives d,double r,double T,double sigma) : model(d,r,T,sigma) { }

    double normalCDF(double value) {
        return 0.5 * erfc(-value / sqrt(2));
    }

    double CalculPrice() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        return (d.get_S0()*normalCDF(d1) - normalCDF(d2)*d.get_K()*exp(-r*T) );
    }

    double CalculDelta() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        return (normalCDF(d1));
    }

    double CalculGamma() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        return (exp(-d1*d1/2)/sqrt(2*M_PI) /(d.get_S0()*sigma*sqrt(T)) ) ;
    }

    double CalculVega() {
        double d1 = (log(d.get_S0() / d.get_K()) + (r + sigma*sigma / 2)*T) / (sigma*sqrt(T));
        cout << d1 << endl ;
        return(exp(-d1*d1/2)/sqrt(2*M_PI) * d.get_S0() *sqrt(T)) ;
    }
};

class MC_Asian_Call : public model {
	int N; // number of path incrememnts, to be input as a model parameter
public :
	MC_Asian_Call(Derivatives d, double r, double T, double sigma, int N) : model(d, r, T, sigma) {
		// constructor to include N in additional to existing parameters
		this -> N = N;
	}

	double geo_average(vector<double> &vec) {
		// function to caluclate geometric average, given vector, uses finds the average of the log values first, then gives exp( average(log values) )
		double log_sum = 0.0; 
		for (unsigned int i = 0; i < vec.size(); i++) {
			log_sum = log_sum + log(vec[i]);
		}
		return exp(log_sum / vec.size());
	}

	vector<double> generate_path() {
		double h = T / N;
		vector<double> path(N+1);
		vector<double> Z = generate_normal(0.0, 1.0, N);

		path[0] = d.get_S0();

		for (int i = 0; i < N; i++) {
			path[i+1] = path[i] + r * path[i]*h + sigma * path[i]*sqrt(h)*Z[i];
		}
		return path;
	}

	double CalculPrice(int M) {
		// does not use Derivative d as input as it is present in member data;
		// M is number of MC simulations;

		// payoff for asian option is: A_T = ( geo_average( price_path ) - K )^+

		vector<double> A_T_vec;
		for (int i = 0; i < M; i++) {
			// M times of MC simulations
			vector<double> price_path = generate_path();
			// cout << "PATH[1]: " << price_path[1] << endl;
			double average = geo_average(price_path);
			// cout << "AVERAGE:" << average << endl;
			// cout << "indicator = " << ((average - d.get_K()) >= 0 )<< endl;
			// cout << "A_T = " << ((average - d.get_K()) >= 0) * (average - d.get_K()) << endl;
			A_T_vec.push_back(((average - d.get_K()) >= 0) * (average - d.get_K()));
		}

		return mean_vector(A_T_vec);
	}

};


class CF_Asian_Call : public model {
public:
	CF_Asian_Call(Derivatives d, double r, double T, double sigma) : model(d, r, T, sigma) { }
	double normalCDF(double value) {
		// PROPOSAL: move this higher up as a public function instead of a member function?
		return 0.5 * erfc(-value / sqrt(2));
	}

	double CalculPrice() {
		
		double t = 0.0; // initial time
		double G_t = d.get_S0();
		double mu_bar = (r - sigma * sigma / 2) * pow(T - t, 2) / (2 * T);
		double sigma_bar = sqrt(sigma*sigma / (T*T) * pow(T - t, 3) / 3);
		double d2 = 1.0 / sigma_bar * (t / T * log(G_t) + (T - t) / T * log(d.get_S0()) + mu_bar - log(d.get_K()));
		double d1 = d2 + sigma_bar;

		return exp(-r * (T - t)) * (pow(G_t, t / T) * pow(d.get_S0(), (T - t) / T) * exp(mu_bar + pow(sigma_bar, 2) / 2) * normalCDF(d1) - d.get_K()*normalCDF(d2));
	}

};

int main() {
    Derivatives call(100,100);
    MonteCarlo mc(call ,0.05,1.0,0.4, 100) ;
    Black_Scholes bs(call ,0.05,1.0,0.4) ;
    // price of the option
    double price_call_mc = mc.CalculPrice(call, 1000000);
    cout << "The Monte_Carlo price of the call is equal to "<< price_call_mc << endl;
    double price_call_bs = bs.CalculPrice();
    cout << "The BS price of the call is equal to "<< price_call_bs << endl;

    //mean
    double variance_MC=mc.CalculVariance(call,1000000);
    cout << "The Monte_Carlo variance of the call is equal to "<< variance_MC << endl;


    // delta
    double delta_mc = mc.CalculDelta(call, 1000000, 0.1);
    cout << "The Monte_Carlo delta is equal to "<< delta_mc << endl;
    double delta_bs = bs.CalculDelta();
    cout << "The BS delta is equal to "<< delta_bs << endl;
    double delta_new=mc.CalculDelta_2(call,1000000);
    cout << "The Monte_Carlo delta_1 is equal to "<< delta_new << endl;
    // gamma
    double gamma_mc = mc.CalculGamma(call, 1000000, 0.001);
    cout << "The MonteCarlo gamma is equal to "<< gamma_mc << endl;
    double gamma_bs = bs.CalculGamma();
    cout << "The BS gamma is equal to "<< gamma_bs << endl;
    // vega
    double vega_bs = bs.CalculVega();
    cout << "The vega is equal to " << vega_bs << endl;


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

	Derivatives asian_call(100.0, 100.0);
	MC_Asian_Call mc_asian_call(asian_call, 0.05, 1.0, 0.4, 1000);
	CF_Asian_Call cf_asian_call(asian_call, 0.05, 1.0, 0.4);

	// price of the option
	cout << endl;
	double asian_call_price_mc = mc_asian_call.CalculPrice(10000);
	double asian_call_price_cf = cf_asian_call.CalculPrice();
	cout << "The Monte Carlo price of the asian call is equal to " << asian_call_price_mc << endl;
	cout << "The Closed form price of the asian call is equal to " << asian_call_price_cf << endl;



	int dummy;
	cin >> dummy;
    return 0;
}
 