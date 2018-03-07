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
#include <limits>

#include "Normal.h"

using namespace std;

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

double mean( vector<double> &vec ) {
    double sum = 0;
    int n = vec.size();
    for( unsigned int i = 0; i < n; i++) {
        sum += vec[i];
    }
    return sum / n;
}

double objective_F( vector<double> &a_vec, double x, double K ) {
    return a_vec[0] + a_vec[1] * x + a_vec[2] * pow(x,2) - ( x - K ) * ( x - K > 0 );
}

double objective_f( vector<double> &a_vec, double x, double K ) {
    return a_vec[1] + 2 * a_vec[2] * x - ( x - K > 0 );
}

double dot_product( vector<double> &x, vector<double> &y ) {
    double sum = 0;
    for( unsigned int i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

void print( vector<vector<double> > &A ) {
    for( unsigned int i = 0; i < A.size(); i++) {
        for( unsigned int
            j = 0; j < A[0].size(); j++ ) {
            cout << A[i][j] << ", ";
        }
        cout << endl;
    }
}

void print( vector<double> &x ) {
    for( unsigned int j = 0; j < x.size(); j++ ) {
        cout << x[j] << ", ";
    }
    cout << endl;
}

vector<vector<double> > transpose( vector<vector<double> > &A ) {
    vector<vector<double> > res;
    vector<double> row;
    
    for( unsigned int j = 0; j < A[0].size(); j++ ) {
        row.clear();
        for( unsigned int i = 0; i < A.size(); i++ ) {
            row.push_back( A[i][j] );
        }
        res.push_back(row);
    }
    
    return res;
}

vector<vector<double> > matrix_product( vector<vector<double> > &A, vector<vector<double> > &B ) {
    
    vector<vector<double> > res, B_T;
    vector<double> row;
    
    B_T = transpose( B );
    
    for( unsigned int i = 0; i < A.size(); i++ ) {
        row.clear();
        for( unsigned int j = 0; j < B[0].size(); j++ ) {
            row.push_back( dot_product( A[i], B_T[j] ) );
        }
        res.push_back( row );
    }

    return res;
}

vector<vector<double> > inverse_3x3( vector<vector<double> > &m ) {
    
    double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    
    if( det == 0 ) {
        cout << "WARNING: singular matrix" << endl;
        print( m );
    }
    
    double invdet = 1 / det;
    
    vector<vector<double> > res;
    vector<double> row;

    row.clear();
    row.push_back( (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet );
    row.push_back( (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet );
    row.push_back( (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet );
    res.push_back(row);
    
    row.clear();
    row.push_back( (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet );
    row.push_back( (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet );
    row.push_back( (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet );
    res.push_back(row);
    
    row.clear();
    row.push_back( (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet );
    row.push_back( (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet );
    row.push_back( (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet );
    res.push_back(row);
    
    return res;
}

vector<vector<double> > X_from_vec( vector<double> &x ) {
    vector<vector<double> > X;
    vector<double> row;
    for( unsigned int i = 0; i < x.size(); i++) {
        row.clear();
        row.push_back(1);
        row.push_back(x[i]);
        row.push_back( pow(x[i],2) );
        
        X.push_back( row );
    }
    return X;
}

vector<vector<double> > xTx_inverse( vector<double> &x ) {
    vector<vector<double> > res;
    
    double M = x.size();
    double sum_1 = 0, sum_2 = 0, sum_3 = 0, sum_4 = 0;
    
    for( unsigned int i = 0; i < M; i++ ) {
        sum_1 += x[i];
        sum_2 += pow(x[i],2);
        sum_3 += pow(x[i],3);
        sum_4 += pow(x[i],4);
    }
    
    vector<double> row;
    row.push_back(M);
    row.push_back(sum_1);
    row.push_back(sum_2);
    res.push_back(row);
    
    row.clear();
    row.push_back(sum_1);
    row.push_back(sum_2);
    row.push_back(sum_3);
    res.push_back(row);
    
    row.clear();
    row.push_back(sum_2);
    row.push_back(sum_3);
    row.push_back(sum_4);
    res.push_back(row);
    
    res = inverse_3x3(res);
    return res;
    
}

vector<vector<double> > xTy ( vector<double> &x, vector<double> &y ) {
    vector<vector<double> > res;
    vector<double> row;
    double sum_1 = 0, sum_2 = 0, sum_3 = 0;
    
    for( unsigned int i = 0; i < x.size(); i++ ) {
        sum_1 += y[i];
        sum_2 += x[i] * y[i];
        sum_3 += pow( x[i], 2 ) * y[i];
    }
    
    row.push_back( sum_1 );
    res.push_back(row);
    
    row.clear();
    row.push_back( sum_2 );
    res.push_back(row);
    
    row.clear();
    row.push_back( sum_3 );
    res.push_back(row);
    
    return res;
}

vector<double> a_estimate( vector<double> &x, vector<double> &y ) {
    vector<vector<double> > XTX, res, XTY;
    
    XTX = xTx_inverse( x );
    XTY = xTy(x,y);
    
    res = matrix_product(XTX,XTY);
    res = transpose(res);
    
    return res[0];
}

double solve_boundary( vector<double> &x, vector<double> &y, double K ) {
    vector<double> a = a_estimate(x,y);
    
    if( pow(a[1]+1,2) -4*a[2]*(a[0]-K) >= 0 ) {
        return ( -(a[1]+1) + pow( pow(a[1]+1,2) -4*a[2]*(a[0]-K), 0.5) ) / ( 2 * a[2] );
    } else {
        // if no intersection => V_C is always larger, so boundary is infinity;
        return 0;
    }
}

//Normal random number generator
Normal normal(Custom);


class american_moving_average_asian {
    double T; // terminal time
    unsigned int N; // number of time partitions
    double K; // strik
    double W; // maximum width of the moving average
    bool is_call;
    
    vector<vector<double> > Z, A;
    vector<double> b;
    
    
    vector<vector<double> > compute_S( double S0, double r, double v ) {
        
        vector<vector<double> > S;
        vector<double> path;
        
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        
        for( unsigned int i = 0; i < Z.size(); i++ ) {
            // initialize the ith path
            path.clear();
            
            // initialize price at time 0 = S0;
            path.push_back( S0 );
            
            // calculate the rest of the path;
            for (unsigned int j = 0; j < Z[0].size(); j++ ) {
                path.push_back( path[i] + r * path[i]*dt + v * path[i]*dt_sqr*Z[i][j] ); // + 0.5 * v * ( v * path[i] ) * dt * ( pow(Z[i][j],2) - 1 ) );
            }
            
            S.push_back( path );
        }
        return S;
    }
    
    void compute_A( double S0, double r, double v ) {
        A.clear();
        
        vector<vector<double> > S = compute_S( S0, r, v );
        vector<double> path;
        double log_sum;
        
        for( unsigned int i = 0; i < S.size(); i++ ) {
            // initialize the ith path of the geometric average
            path.clear();
            
            for(unsigned int j = 0; j < W; j++ ) {
                // the first W geometric average dont have the full window
                log_sum = 0;
                for( unsigned l = 0; l < j+1; l++ ) {
                    // the jth < W geometric averages have average of j+1 terms
                    log_sum += log( S[i][l] );
                }
                path.push_back( exp( log_sum/(j+1) ) );
            }
            
            for(unsigned int j=W; j < S[0].size(); j++ ) {
                // the remaining geometric averages have full winodw W
                log_sum = 0;
                
                for( unsigned l = j-W+1; l <= j; l++ ) {
                    // the jth geometric has first term at j-W+1, last term at j, adding up to W terms
                    log_sum += log( S[i][l]) ;
                }
                path.push_back( exp( log_sum/(W) ) );
            }
            A.push_back( path );
            
        }
    }
    
    void compute_Z( int no_sims ) {
        Z.clear();
        vector<double> row;
        for(unsigned int i = 0; i< no_sims; i++ ) {
            row = normal.generate(N);
            Z.push_back( row );
        }
    }
    
    double V( int m_id, int t_id, double &r ) {
        if( t_id == N ) {
            // This is the last payoff in the path, there is no V_c to compare with
            return ( K - A[m_id][t_id] ) * (  K - A[m_id][t_id] > 0 );
        } else if ( t_id < N ) {
            if( A[m_id][t_id] <= b[t_id] ) {
                return K - A[m_id][t_id];
            } else {
                return exp(-r*T) * V( m_id, t_id + 1, r);
            }
        } else {
            cout << "ERROR @compute Z" << endl;
            return -1;
        }
    }
    
    void compute_b(double r) {
        b.clear();
        b.resize( A[0].size() - 1 );
        vector<double> y, x;
        
        for(int t_id = (N-1); t_id >= (1); t_id--) {
            
            // clean up vectors x,y
            y.clear();
            x.clear();
            
            for( unsigned int i = 0; i < A.size(); i++ ) {
                y.push_back( V(i, t_id + 1, r) * exp(-r*T) );
                x.push_back( A[i][t_id] );
            }
            
            b[t_id] = solve_boundary( x, y, K );
        }
    }
    
    
    void MC_euler_pricing(double S0, double r, double v, unsigned int no_sims, vector<double>& res) {
        // prameter initlisation
        clock_t c;
        c = clock();
        double duration;
        
        compute_Z( no_sims );
        compute_A( S0, r, v); // if need to change parameter (i.e. finite difference greeks), can change at this step
        
        compute_b(r);
        
        // the payoff of immediate exercies at time 0
        double V_E = ( K - S0 ) * ( K - S0 >= 0 );
        
        // summing the payoff for holding the put option
        double sum = 0;
        for( unsigned int i = 0; i < no_sims; i++ ) {
            sum += V(i, 1, r);
            cout << V(i, 1, r) << endl;
        }
        double V_C = exp(-r*T) * sum / no_sims;
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        
        print(b);
        
        // return results
        res.push_back( max( V_E, V_C ) ); // price
        res.push_back( duration ); //time
    }
    
    
public:
    
    // an asian option has time to maturity T,K,initial price and interest rate
    american_moving_average_asian(const string& type,double T,int N, double K, int W):T(T),N(N),K(K),W(W){
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
            MC_euler_pricing(S0, r, v, no_sims, res);
        } else if (icompare(method,"milstein")){
            // MC_milstein_pricing(S0,r,v,no_sims, res);
        } else if (icompare(method,"analytic")){
            // analytic_solution_pricing(S0, r, v, res);
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
    
};


int main() {
    vector<double> a_vec, s_vec, x_vec, y_vec;
    string type = "call";
    int T = 1;
    unsigned int N = 100;
    double K = 100;
    double W = 2;
    double S0 = 100;
    double r = 0.05;
    double v = 0.4;
    int no_sims = 1000;
    
    srand(time(NULL));
    american_moving_average_asian opt(type,T,N,K,W);
    
    opt.calculate_price("euler", S0, r, v, no_sims, true);
    
}
