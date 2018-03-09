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
    
    double sum_1 = 0, sum_2 = 0, sum_3 = 0, sum_4 = 0;
    
    for( unsigned int i = 0; i < x.size(); i++ ) {
        sum_1 += x[i];
        sum_2 += pow(x[i],2);
        sum_3 += pow(x[i],3);
        sum_4 += pow(x[i],4);
        
        /*
        cout << i << endl;
        cout << x[i]        << " added to " << sum_1 << endl;
        cout << pow(x[i],2) << " added to " << sum_2 << endl;
        cout << pow(x[i],3) << " added to " << sum_3 << endl;
        cout << pow(x[i],4) << " added to " << sum_4 << endl;
        cout << endl;
        */
    }
    
    vector<double> row;
    row.push_back(x.size());
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
    
    
    
    // cout << "   XTX = " << endl;
    // print( res );
    
    res = inverse_3x3(res);
    return res;
    
}

vector<vector<double> > xTy ( vector<double> &x, vector<double> &y ) {
    vector<vector<double> > res;
    vector<double> row;
    double sum_1 = 0, sum_2 = 0, sum_3 = 0;
    
    for( unsigned int i = 0; i < x.size(); i++ ) {
        sum_1 += y[i];
        sum_2 += y[i] * x[i];
        sum_3 += y[i] * pow( x[i], 2 );
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
    // cout << "   (XTX)inverse = " << endl;
    // print(XTX);
    
    XTY = xTy(x,y);
    // cout << "   XTY = " << endl;
    // print(XTY);
    
    res = matrix_product(XTX,XTY);
    res = transpose(res);
    
    return res[0];
}

vector<double> delta_F( vector<double> &x, vector<double> &y, vector<double> &a ) {
    vector<double> res(3,0);
    
    for( unsigned int j = 0; j <x.size(); j++ ) {
        res[0] += 2 * ( y[j] - a[0] - a[1] * x[j] - a[2] * pow(x[j],2) ) * ( -1 );
        res[1] += 2 * ( y[j] - a[0] - a[1] * x[j] - a[2] * pow(x[j],2) ) * ( -x[j] );
        res[2] += 2 * ( y[j] - a[0] - a[1] * x[j] - a[2] * pow(x[j],2) ) * ( -x[j]*x[j] );
    }
    return res;
}

double F( vector<double> &x, vector<double> &y, vector<double> &a ) { \
    double res = 0;
    
    for( unsigned int j = 0; j <x.size(); j++ ) {
        res += pow(( y[j] - a[0] - a[1] * x[j] - a[2] * pow(x[j],2) ), 2);
    }
    return res;
}

double gd_gamma( vector<double> a_old, vector<double> a_new, vector<double> delta_old, vector<double> delta_new ) {
    double numerator = 0, denomiator = 0;
    
    for( unsigned int i = 0; i < a_old.size(); i++ ) {
        numerator += ( a_new[i] - a_old[i] ) * ( delta_new[i] - delta_old[i] );
        denomiator += pow( delta_new[i] - delta_old[i], 2 );
    }
    //cout << numerator << endl;
    //cout << denomiator << endl;
    return abs( numerator / denomiator );
}

double norm( vector<double> x ) {
    double sum = 0;
    for( unsigned int i = 0; i <x.size(); i++ ) {
        sum += pow(x[i],2);
    }
    return pow(sum,0.5);
}

vector<double> a_gd_estimate( vector<double> &x, vector<double> &y ) {
    vector<double> a_old(3,100), a_new(3,0), delta_F_new(3,0), delta_F_old(3,0);
    double gamma = 1, F_val_old, F_val_new;
    bool done = false;
    double no_iter = 1000000;
    double tol = pow(10,-2);
    double diff = 1;
    int i = 0;
    
    hline(50);
    
    cout << "x = ";
    print(x);
    cout << "y = ";
    print(y);
    cout << endl;
    
    while( diff > tol ) {
        // apply gradient desecnt with an initial gamma = 1
        delta_F_new = delta_F( x,y, a_old );
        F_val_old = F(x,y,a_old);
        a_new[0] = a_old[0] - gamma * delta_F_new[0];
        a_new[1] = a_old[1] - gamma * delta_F_new[1];
        a_new[2] = a_old[2] - gamma * delta_F_new[2];
        F_val_new = F(x,y,a_new);
        
        while( F_val_old < F_val_new ) {
            // decrecase gamma by a factor of 10 if new F is not smaller than old F
            // cout << "DECREASING gamma" << endl;
            gamma = gamma * 0.5;
            a_new[0] = a_old[0] - gamma * delta_F_new[0];
            a_new[1] = a_old[1] - gamma * delta_F_new[1];
            a_new[2] = a_old[2] - gamma * delta_F_new[2];
            F_val_new = F(x,y,a_new);
        }
        
        
        //cout << "i = " << i << endl;
        /*
         cout << "   a_old = ";
         print(a_old);
         cout << "   F_old = " << F_val_old << endl;
         cout << "   delta F = ";
         print( delta_F_new );
         cout << "   a_new = ";
         print(a_new);
         cout << "   F_new = " << F_val_new << endl;
         */
        //cout << "   gamma used = " << gamma << endl;
        //cout << "   delta norm = " << norm(delta_F_new) << endl;
        //cout << endl;
        
        diff = norm(delta_F_new);
        // renew variables for next step
        delta_F_old = delta_F_new;
        a_old = a_new;
        // gamma = abs( gd_gamma(a_old, a_new, delta_F_old, delta_F_new));
        gamma = 1;
        i += 1;
        
        if( i > no_iter ) {
            break;
        }
    }
    if( diff < tol) {
        cout << "Required tolerance of delta norm " << tol << " reached at i = " << i << endl;
        cout << endl;
        done = true;
        // break;
    }
    
    
    if(!done) {
        cout << "Required tolerance has not been reached at i = " << no_iter << " still at " << norm(delta_F_new) << endl;
    }
    
    
    return a_new;
}

double solve_boundary( vector<double> &x, vector<double> &y, double K ) {
    vector<double> a = a_estimate(x,y); // a_gd = a_gd_estimate(x,y);
    
    // print(a);
    
    if( pow(a[1]+1,2) -4*a[2]*(a[0]-K) >= 0 ) {
        return ( -(a[1]+1) + pow( pow(a[1]+1,2) -4*a[2]*(a[0]-K), 0.5) ) / ( 2 * a[2] );
    } else {
        // if no intersection => V_C is always larger, so boundary is 0;
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
    
    
    vector<vector<double> > compute_S( string method, double S0, double r, double v ) {
        
        vector<vector<double> > S;
        vector<double> path;
        
        double dt = (double)T / N;
        double dt_sqr = pow((double)T / N, 0.5); // ADDED: Casting T to double
        
        if( icompare(method, "euler") ) {
            for( unsigned int i = 0; i < Z.size(); i++ ) {
                // initialize the ith path
                path.clear();
                
                // initialize price at time 0 = S0;
                path.push_back( S0 );
                
                // calculate the rest of the path;
                for (unsigned int j = 0; j < Z[0].size(); j++ ) {
                    path.push_back( path[j] + r * path[j]*dt + v * path[j]*dt_sqr*Z[i][j] );
                }
                
                S.push_back( path );
            }
        } else {
            cout << "WARNING: unclear method of " << method << endl;
        }
        
        return S;
    }
    
    void compute_A( string method,  double S0, double r, double v ) {
        A.clear();
        
        vector<vector<double> > S = compute_S( method, S0, r, v );
        /*
        cout << "S = " << endl;
        print(S);
         */
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
        // implement the V as definined in American put option lecturue notes
        
        double multiplier = 1;
        double dt = (double) T/N;
        
        // cout << "FIRST STEP = " << t_id << endl;
        
        while( t_id < N && A[m_id][t_id] > b[t_id] ) {
            // 1st condition:   if t_id < N mean it is not the last time step
            // 2nd condition:   if A[m_id][t_id] > b[t_id], means not exercise ad this time step
            //                  then increase one time step, and also mulitiple a discount
            t_id += 1;
            multiplier = multiplier * exp(-r*dt);
        }
        
        if( t_id == N ) {
            // This is the last payoff in the path, there is no V_c to compare with
            // cout << "   LAST STEP = " << t_id << endl;
            return multiplier * ( K - A[m_id][t_id] ) * (  K - A[m_id][t_id] > 0 );
        } else {
            
            // cout << "   LAST STEP = " << t_id << endl;
            return multiplier * ( K - A[m_id][t_id] );
        }
    }
    
    
    vector<double> V_detailed( int m_id, int t_id, double &r ) {
        // mostly a copy of function V, but output more stats;
        
        double multiplier = 1;
        double dt = (double) T/N;
        
        int FIRST = t_id;
        // cout << "FIRST STEP = " << t_id << endl;
        
        while( t_id < N && A[m_id][t_id] > b[t_id] ) {
            // 1st condition:   if t_id < N mean it is not the last time step
            // 2nd condition:   if A[m_id][t_id] > b[t_id], means not exercise ad this time step
            //                  then increase one time step, and also mulitiple a discount
            t_id += 1;
            multiplier = multiplier * exp(-r*dt);
        }
        
        int LAST = t_id;
        vector<double> res;
        
        if( t_id == N ) {
            // This is the last payoff in the path, there is no V_c to compare with
            // cout << "   LAST STEP = " << t_id << endl;
            res.push_back( multiplier * ( K - A[m_id][t_id] ) * (  K - A[m_id][t_id] > 0 ) );
            res.push_back( FIRST );
            res.push_back( LAST );
            return res;
        } else {
            // cout << "   LAST STEP = " << t_id << endl;
            res.push_back( multiplier * ( K - A[m_id][t_id] ) ) ;
            res.push_back( FIRST );
            res.push_back( LAST );
            return res;
        }
    }
    
    void compute_b(double r) {
        b.clear();
        b.resize( A[0].size() - 1 );
        vector<double> y, x;
        
        double dt = (double) T/N;
        for(int t_id = (N-1); t_id >= (1); t_id--) {
            
            // clean up vectors x,y
            y.clear();
            x.clear();
            
            for( unsigned int i = 0; i < A.size(); i++ ) {
                y.push_back( V(i, t_id + 1, r) * exp(-r*dt) );
                x.push_back( A[i][t_id] );
            }
            
            // cout << "t_id = " << t_id << endl;
            double b_val = solve_boundary( x, y, K );
            if( !isnan(b_val) ) {
                b[t_id] = b_val;
            } else {
                cout << "ERROR: b has NAN at t_id = " << t_id << endl;
                break;
            }
            // cout << "t_id = " << t_id << endl;
            // cout << "   b = ";
            // print(b);
            // cout << endl;
        }
    }
    
    
    void MC_pricing(string method, double S0, double r, double v, unsigned int no_sims, vector<double>& res) {
        // prameter initlisation
        clock_t c;
        c = clock();
        double duration;
        double dt = (double) T/N;
        
        compute_Z( no_sims );
        // cout << "Z = " << endl;
        // print(Z);
        
        compute_A( method, S0, r, v); // if need to change parameter (i.e. finite difference greeks), can change at this step
        // cout << "A = " << endl;
        // print(A);
        
        compute_b(r);
        // cout << "b = ";
        // print(b);
        
        
        // the payoff of immediate exercies at time 0
        double V_E = ( K - S0 ) * ( K - S0 >= 0 );
        
        // summing the payoff for holding the put option
        double sum = 0;
        vector<double> path_res;
        for( unsigned int i = 0; i < no_sims; i++ ) {
            path_res = V_detailed(i, 1, r);
            sum += path_res[0];
            if( path_res[2] < N ) {
                // cout << "Path " << i+1 << " early exercise at " << path_res[2] << " with price " << abs(path_res[0]) * exp(-r*dt) << endl;
            } else {
                // cout << "Path " << i+1 << " exercise at last step with price " << abs(path_res[0]) * exp(-r*dt) << endl;
            }
        }
        double V_C = exp(-r*dt) * sum / no_sims;
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        
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
        
        //Use selected method to calcuate c_0 and put it in res;
        MC_pricing( method, S0, r, v, no_sims, res );
        
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
    unsigned int N = 10;
    double K = 100;
    double W = 2;
    double S0 = 100;
    double r = 0.05;
    double v = 0.4;
    double no_sims = 100;
    
    srand(time(NULL));
    american_moving_average_asian opt(type,T,N,K,W);
    for( unsigned int i = 0; i < 5; i++ ) {
        opt.calculate_price("euler", S0, r, v, no_sims*pow(10,i), true);
    }
}
;
