#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double mean( vector<double> &vec ) {
    double sum = 0;
    int n = vec.size();
    for( unsigned int i = 0; i < n; i++) {
        sum += vec[i];
    }
    return sum / n;
}

double objective_F( vector<double> &a_vec, vector<double> s_vec, double x, double K ) {
    s_vec.push_back(x);
    return a_vec[0] + a_vec[1] * x + a_vec[2] * pow(x,2) - ( mean(s_vec) - K ) * ( mean(s_vec) - K > 0 );
}

double objective_f( vector<double> &a_vec, vector<double> s_vec, double x, double K ) {
    int w = s_vec.size() + 1;
    s_vec.push_back(x);
    return a_vec[1] + 2 * a_vec[2] * x - ( x / w ) * ( mean(s_vec) - K > 0 );
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

vector<vector<double> > xTx_inverse( vector<vector<double> > &X ) {
    vector<vector<double> > res, XT;
    
    XT = transpose(X);
    res = matrix_product( XT, X );
    res = inverse_3x3(res);
    return res;
    
}

vector<double> a_estimate( vector<double> &x, vector<double> &y ) {
    vector<vector<double> > XTX, X, XT, res, Y;
    
    X = X_from_vec( x );
    XT = transpose(X);
    XTX = xTx_inverse( X );
    Y.push_back(y);
    Y = transpose(Y);
    
    res = matrix_product(XTX,XT);
    res = matrix_product(res,Y);
    res = transpose(res);
    
    return res[0];
}

double solve_boundary( vector<double> &x, vector<double> &y, vector<double> s_vec, double K ) {
    
    vector<double> a_vec = a_estimate(x,y);
    
    double req_tol = 0.001, diff = 1;
    
    int counter = 1;
    
    double new_value = s_vec.back() ; // use last price in the moving window as initial guess value
    s_vec.pop_back();
    double old_value;
    
    while( abs(diff) > req_tol ) {
        
        old_value = new_value;
        cout << "Counter: " << counter++ << ", ";
        
        diff = objective_F( a_vec, s_vec, old_value, K ) / objective_f( a_vec, s_vec, old_value, K);
        cout << "diff = " << diff << ", " ;
        
        new_value = old_value - diff;
        
        cout << new_value << endl;
    }
    cout << "DONE! with diff = " << diff << endl;
    return new_value;
}

int main() {
    vector<double> a_vec, s_vec, x_vec, y_vec;
    double K = 0;

    
    x_vec.push_back(2);
    x_vec.push_back(3);
    x_vec.push_back(0);
    
    y_vec.push_back(2);
    y_vec.push_back(3);
    y_vec.push_back(0);
    
    a_vec = a_estimate(x_vec,y_vec);
    print(a_vec);
    
    s_vec.push_back(10);
    
    double b = solve_boundary( x_vec, y_vec, s_vec, 0 );
    cout << endl;
    cout << b;
    
}
