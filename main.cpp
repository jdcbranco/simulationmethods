#include <iostream>
#include "Normal.h"
#include "Vanilla.h"
#include "MCModel.h"
#include "BSVanilla.h"
#include "Asian.h"
#include "BSAsian.h"

using namespace std;

Normal normal(Standard);

int main() {
    double s0 = 100.0;
    double strike = 100;
    double sigma = 0.4;
    double r = 0.05;

    vector<int> number_simulations = {1000,5000,10000,25000,50000,100000,200000,500000} ;
    //vector<int> number_simulations = {1000,5000};

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
        ModelResult mcModelResult = mcModel.simulate(simulator,vanilla_call_fd,i,mcModel.getSolver() == Explicit ? 1 : 5);
        //ModelResult mcModelResult = mcModel.simulate(simulator,i,mcModel.getSolver() == Explicit ? 1 : 5, SensitivityMethod::Greeks_by_LR);
        cout << "-------------" << endl;
        cout << "Simulations: " << i << endl;
        cout << mcModelResult;
    }
    //Asian call
    double path_size = 100;
    AsianCall asianCall(strike, 1.0);
    BSAsianCallModel bsAsianModel(asianCall, s0, sigma, r, path_size);
    ModelResult bsAsianModelResult = bsAsianModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes Asian Call: " << endl;
    cout << bsAsianModelResult;
    MCModel<AsianCall> asianMcModel(asianCall, s0, sigma, r, 0.005, Explicit); //Optionally, can try Euler as well. Both work fine.
    Greeks_by_FD::CentralDifferencesSensitivityModel<AsianCall> asian_call_fd(asianMcModel, 0.005);
    Greeks_by_PD::AsianCallSensitivityModel asian_call_pd(asianMcModel);
    Greeks_by_LR::AsianSensitivityModel<AsianCall> asian_call_lr(asianMcModel);
    ModelResult asianMcModelResultFD = asianMcModel.simulate(simulator, asian_call_fd, 100000, path_size);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and "<< path_size <<" steps (Finite Differences): " << endl;
    cout << asianMcModelResultFD;
    ModelResult asianMcModelResultPD = asianMcModel.simulate(simulator, asian_call_pd, 100000, path_size);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and "<< path_size <<" steps (Pathwise Differentiation): " << endl;
    cout << asianMcModelResultPD;
    ModelResult asianMcModelResultLR = asianMcModel.simulate(simulator, asian_call_lr, 100000, path_size);
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