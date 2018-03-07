#include <iostream>
#include "Normal.h"
#include "Vanilla.h"
#include "MCModel.h"
#include "BSVanilla.h"
#include "Asian.h"

using namespace std;

Normal normal(Standard);

int main() {
    double s0 = 100.0;
    double strike = 100;
    double sigma = 0.4;
    double r = 0.05;

    vector<int> number_simulations = {1000,5000,10000,25000,50000,100000,200000,500000} ;

    //European call
    VanillaCall vanillaCall(strike, 1.0);
    BSCallModel bsModel(vanillaCall, s0, sigma, r);
    MCModel mcModel(vanillaCall, s0, sigma, r, 0.005, Explicit);
    Simulator simulator(normal,true); //True=Antithetic
    ModelResult bsModelResult = bsModel.calculate(); //This is tested and matches the existing result
    cout << "-------------" << endl;
    cout << "Black Scholes: " << endl;
    cout << bsModelResult;
    for(int i: number_simulations) {
        ModelResult mcModelResult = mcModel.simulate(simulator,i,mcModel.getSolver() == Explicit ? 1 : 5);
        cout << "-------------" << endl;
        cout << "Simulations: " << i << endl;
        cout << mcModelResult;
    }
    //Asian call
    GFixedStrikeAsianCall asianCall(strike, 1.0);
    MCModel asianMcModel(asianCall, s0, sigma, r, 0.005, Explicit); //Optionally, can try Euler as well. Both work fine.
    ModelResult asianMcModelResult = asianMcModel.simulate(simulator,100000, 100);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and 100 steps: " << endl;
    cout << asianMcModelResult;

    return 0;
}