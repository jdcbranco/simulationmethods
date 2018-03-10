#include <iostream>
#include <map>
#include <fstream>

#include "Normal.h"
#include "Vanilla.h"
#include "MCModel.h"
#include "BSVanilla.h"
#include "Asian.h"
#include "BSAsian.h"

using namespace std;

Normal normal(Standard);

pair<double,double> price_mean_variance(vector<ModelResult> result){
    double sum = 0.0;
    double sum_squares = 0.0;
    double size = result.size();
    for (int i = 0; i < result.size(); ++i) {
        sum += result[i].getPrice();
    }
    double mean = sum / size;
    for (int i = 0; i < result.size(); ++i) {
        sum_squares += pow(result[i].getPrice()- mean,2.0);
    }

    double variance = sum_squares / (size-1);
    return pair<double,double>(mean,variance);
}

enum class Greek { Delta, Gamma, Vega };

pair<double,double> greeks_mean_variance(vector<ModelResult> result, Greek greek){
    double sum = 0.0;
    double sum_squares = 0.0;
    double size = result.size();
    for (int i = 0; i < result.size(); ++i) {
        switch (greek) {
            case Greek::Delta:
                sum += result[i].getDelta();
                break;
            case Greek::Gamma:
                sum += result[i].getGamma();
                break;
            case Greek::Vega:
                sum += result[i].getVega();
                break;
        }
    }
    double mean = sum / size;
    for (int i = 0; i < result.size(); ++i) {
        switch (greek) {
            case Greek::Delta:
                sum_squares += pow(result[i].getDelta()- mean,2.0);
                break;
            case Greek::Gamma:
                sum_squares += pow(result[i].getGamma()- mean,2.0);
                break;
            case Greek::Vega:
                sum_squares += pow(result[i].getVega()- mean,2.0);
                break;
        }
    }

    double variance = sum_squares / (size-1);
    return pair<double,double>(mean,variance);
}

int main() {
    double s0 = 100.0;
    double strike = 100;
    double sigma = 0.4;
    double r = 0.05;

    vector<int> number_simulations = {1000,5000,10000,25000,50000};//,100000,200000, 300000, 400000 ,500000} ;
    map<int,vector<ModelResult>> results;
    int number_runs = 10;

    //European call
    VanillaCall vanillaCall(strike, 1.0);
    BSCallModel bsModel(vanillaCall, s0, sigma, r);
    MCModel<VanillaCall> mcModel(vanillaCall, s0, sigma, r, 0.005, Explicit);
    Greeks_by_FD::CentralDifferencesSensitivityModel<VanillaCall> vanilla_call_fd(mcModel,0.005);
    Greeks_by_PD::VanillaCallSensitivityModel vanilla_call_pd(mcModel);
    Greeks_by_LR::VanillaSensitivityModel<VanillaCall> vanilla_call_lr(mcModel);
    Simulator simulator(normal,true); //True=Antithetic
    ModelResult bsModelResult = bsModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes European Call: " << endl;
    cout << bsModelResult;
    for(int i: number_simulations) {
        //continue;
        vector<ModelResult> runs;
        cout << "-------------" << endl;
        cout << "Simulations: " << i << endl;
        int nr = i>=100000? 10: number_runs;
        for(int j = 0; j<nr; j++) {
            ModelResult mcModelResult = mcModel.simulate(simulator,vanilla_call_lr,i,mcModel.getSolver() == Explicit ? 1 : 5);
            runs.push_back(mcModelResult);
            cout << "(" << j << ") " << mcModelResult;
        }
        results[i] = runs;
    }
    cout << "-------------" << endl;
    cout << "Prices (Mean and Stdev) for "<<number_runs<< " repeated runs" << endl;
    for(int i: number_simulations) {
        auto statistics = price_mean_variance(results[i]);
        cout << "(" <<i<<") " << statistics.first << " " << sqrt(statistics.second) << endl;
    }
    cout << "-------------" << endl;
    cout << "Delta (Mean and Stdev) for "<<number_runs<<" repeated runs" << endl;
    for(int i: number_simulations) {
        auto statistics = greeks_mean_variance(results[i], Greek::Delta);
        cout << "(" <<i<<") " << statistics.first << " " << sqrt(statistics.second) << endl;
    }
    cout << "-------------" << endl;
    cout << "Gamma (Mean and Stdev) for "<<number_runs<<" repeated runs" << endl;
    for(int i: number_simulations) {
        auto statistics = greeks_mean_variance(results[i], Greek::Gamma);
        cout << "(" <<i<<") " << statistics.first << " " << sqrt(statistics.second) << endl;
    }
    cout << "-------------" << endl;
    cout << "Vega (Mean and Stdev) for "<<number_runs<<" repeated runs" << endl;
    for(int i: number_simulations) {
        auto statistics = greeks_mean_variance(results[i], Greek::Vega);
        cout << "(" <<i<<") " << statistics.first << " " << sqrt(statistics.second) << endl;
    }
    //Asian call
    double asian_dim = 100; //Number of prices considered in the average
    AsianCall asianCall(strike, 1.0, asian_dim);
    ArithmeticAsianCall arithmeticAsianCall(strike, 1.0, asian_dim);
    BSAsianCallModel bsAsianModel(asianCall, s0, sigma, r, asian_dim);
    ModelResult bsAsianModelResult = bsAsianModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes Asian Call: " << endl;
    cout << bsAsianModelResult;
    {
        double path_size = asian_dim;
        cout << "-------------" << endl;
        cout << "Using Explicit Price SDE Solution: " << endl;
        MCModel<ArithmeticAsianCall> asianMcModel(arithmeticAsianCall, s0, sigma, r, 0.005, Explicit);
        Greeks_by_FD::CentralDifferencesSensitivityModel<ArithmeticAsianCall> asian_call_fd(asianMcModel, 0.005);
        ModelResult asianMcModelResultFD = asianMcModel.simulate(simulator, asian_call_fd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Arithmetic Asian Call with 100k paths and " << path_size << " steps (FD): " << endl;
        cout << asianMcModelResultFD;
        //Using Control variate
        asianMcModel.define_control_variate([&](const Path &path) { return asianCall.payoff(path,None); }, bsAsianModelResult.getPrice());
        ModelResult asianMcModelResultWithControlVariate = asianMcModel.simulate(simulator,asian_call_fd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Arithmetic Asian Call with 100k paths and " << path_size << " steps and *Control Variates* (FD): " << endl;
        cout << asianMcModelResultWithControlVariate;
    }
    {
        double path_size = 1;
        cout << "-------------" << endl;
        cout << "Using Explicit Geometric SDE Solution: " << endl;
        MCModel<AsianCall> asianMcModel(asianCall, s0, sigma, r, 0.005, ExplicitGeometricAverage);
        Greeks_by_FD::CentralDifferencesSensitivityModel<AsianCall> asian_call_fd(asianMcModel, 0.005);
        Greeks_by_PD::AsianCallSensitivityModel asian_call_pd(asianMcModel);
        Greeks_by_LR::AsianSensitivityModel<AsianCall> asian_call_lr(asianMcModel);
        ModelResult asianMcModelResultFD = asianMcModel.simulate(simulator, asian_call_fd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Geometric Asian Call with 100k paths and " << path_size << " steps (Finite Differences): " << endl;
        cout << asianMcModelResultFD;
        ModelResult asianMcModelResultPD = asianMcModel.simulate(simulator, asian_call_pd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Geometric Asian Call with 100k paths and " << path_size << " steps (Pathwise Differentiation): " << endl;
        cout << asianMcModelResultPD;
        ModelResult asianMcModelResultLR = asianMcModel.simulate(simulator, asian_call_lr, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Geometric Asian Call with 100k paths and " << path_size << " steps (Likelihood Ratio): " << endl;
        cout << asianMcModelResultLR;

        path_size = asian_dim;
        cout << "-------------" << endl;
        cout << "Using Explicit Price SDE Solution: " << endl;
        MCModel<AsianCall> asianMcModel2(asianCall, s0, sigma, r, 0.005, Explicit);

        double mean_a = (r-sigma*sigma/2)*asianCall.getT()*(asian_dim+1)/(2*asian_dim);
        double var_a = asianCall.getT()*sigma*sigma*((asian_dim+1)*(2*asian_dim+1)/(6.0*asian_dim*asian_dim));
        double control_variate_mean = exp(-r*asianCall.getT())*(s0*exp(mean_a+var_a/2.0));
        asianMcModel2.define_control_variate([&](const Path &path) { return path.geometric_average(None); }, control_variate_mean, 0.8);
        ModelResult asianMcModelResultFDWithControlVariate = asianMcModel2.simulate(simulator, asian_call_fd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Geometric Asian Call with 100k paths and " << path_size << " steps with *Control Variates* (Finite Differences): " << endl;
        cout << asianMcModelResultFDWithControlVariate;

        path_size = 1;
        cout << "-------------" << endl;
        cout << "Using Explicit Geometric SDE Solution: " << endl;

        asianMcModel.define_control_variate([&](const Path &path) { return path.back(None); }, control_variate_mean, 0.7);
        ModelResult asianMcModelResultFDWithControlVariate2 = asianMcModel.simulate(simulator, asian_call_fd, 100000, path_size);
        cout << "-------------" << endl;
        cout << "Geometric Asian Call with 100k paths and " << path_size << " steps with *Control Variates* (Finite Differences): " << endl;
        cout << asianMcModelResultFDWithControlVariate2;
    }
    
    // Normal generator  
    int n_simulations_normal = 10000 ; // number of simulations 
    Normal normal_simulations = Normal(Custom, 0.0,1.0);
    vector<double> normal_simulation = normal_simulations.generate(n_simulations_normal);

    ofstream fout("/output_normal.txt");
    if (! fout.is_open()) { // test that file is open
        cout << "Error opening output file." << endl;
        return -1;
    }

    for (int i = 0 ; i<normal_simulation.size(); i=i+1) {
        fout << normal_simulation[i] << " ";
    }
    fout.close();

    return 0;
}
