#include <iostream>
#include <map>

#include "Normal.h"
#include "Vanilla.h"
#include "MCModel.h"
#include "BSVanilla.h"
#include "Asian.h"
#include "BSAsian.h"

using namespace std;

Normal normal(Custom);

pair<double,double> mean_variance(vector<ModelResult> result){
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

int main() {
    double s0 = 100.0;
    double strike = 100;
    double sigma = 0.4;
    double r = 0.05;

    vector<int> number_simulations = {1000,5000,10000,25000,50000};//,100000,200000, 300000, 400000 ,500000} ;
    map<int,vector<ModelResult>> results;
    int number_runs = 100;

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
        continue;
        vector<ModelResult> runs;
        //ModelResult mcModelResult = mcModel.simulate(simulator,i,mcModel.getSolver() == Explicit ? 1 : 5, SensitivityMethod::Greeks_by_LR);
        cout << "-------------" << endl;
        cout << "Simulations: " << i << endl;
        int nr = i>=100000? 10: number_runs;
        for(int j = 0; j<nr; j++) {
            ModelResult mcModelResult = mcModel.simulate(simulator,vanilla_call_fd,i,mcModel.getSolver() == Explicit ? 1 : 5);
            runs.push_back(mcModelResult);
            cout << "(" << j << ") " << mcModelResult;
        }
        results[i] = runs;
    }
    cout << "-------------" << endl;
    cout << "Prices (Mean and Stdev) for 100 repeated runs" << endl;
    for(int i: number_simulations) {
        auto statistics = mean_variance(results[i]);
        cout << "(" <<i<<") " << statistics.first << " " << sqrt(statistics.second) << endl;
    }
    //Asian call
    double path_size = 1, asian_dim = 100;
    AsianCall asianCall(strike, 1.0, asian_dim);
    BSAsianCallModel bsAsianModel(asianCall, s0, sigma, r, asian_dim);
    ModelResult bsAsianModelResult = bsAsianModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes Asian Call: " << endl;
    cout << bsAsianModelResult;
    MCModel<AsianCall> asianMcModel(asianCall, s0, sigma, r, 0.005, ExplicitGeometricAverage); //Optionally, can try Euler as well. Both work fine.
    Greeks_by_FD::CentralDifferencesSensitivityModel<AsianCall> asian_call_fd(asianMcModel, 0.005);
    Greeks_by_PD::AsianCallSensitivityModel asian_call_pd(asianMcModel);
    Greeks_by_LR::AsianSensitivityModel<AsianCall> asian_call_lr(asianMcModel);
    ModelResult asianMcModelResultFD = asianMcModel.simulate(simulator, asian_call_fd, 1000000, path_size);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and "<< path_size <<" steps (Finite Differences): " << endl;
    cout << asianMcModelResultFD;
    ModelResult asianMcModelResultPD = asianMcModel.simulate(simulator, asian_call_pd, 1000000, path_size);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and "<< path_size <<" steps (Pathwise Differentiation): " << endl;
    cout << asianMcModelResultPD;
    ModelResult asianMcModelResultLR = asianMcModel.simulate(simulator, asian_call_lr, 1000000, path_size);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and "<< path_size <<" steps (Likelihood Ratio): " << endl;
    cout << asianMcModelResultLR;
    //using control variate
//    asianMcModel.define_control_variate([&](const Path &path) { return vanillaCall.payoff(path,None); }, bsModelResult.getPrice());
//    ModelResult asianMcModelResultWithControlVariate = asianMcModel.simulate(simulator,100000, 10);
//    cout << "-------------" << endl;
//    cout << "Asian Call with 100k paths and 10 steps and Control Variates: " << endl;
//    cout << asianMcModelResultWithControlVariate;

    return 0;
}