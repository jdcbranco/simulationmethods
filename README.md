# simulationmethods
Imperial College London Simulation Methods 

To open the project, use CLion (https://www.jetbrains.com/clion/), which is freely available for students.

Building everything should be simple, but if you want, you can run the following cmake command:

"C:\Program Files\JetBrains\CLion 2017.3.1\bin\cmake\bin\cmake.exe" --build C:\Users\jdcbr\CLionProjects\simulationmethods\cmake-build-debug --target all -- -j 2

You should be able to use any cmake binary you may have available.

## Example use of the code

	double sigma = 0.4;
	double r = 0.05;

	vector<int> number_simulations = {1000,5000,10000,25000,50000,100000,200000,500000} ;

    //New code for European call
    VanillaCall vanillaCall(strike, 1.0);
    BSCallModel bsModel(vanillaCall, 100.0, sigma, r);
    MCModel mcModel(vanillaCall, 100.0, sigma, r, 0.005, Explicit);
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
    //New code for Asian call
    AsianCall asianCall(strike, 1.0);
    MCModel asianMcModel(asianCall, 100.0, sigma, r, 0.005, Explicit); //Optionally, can try Euler as well. Both work fine.
    ModelResult asianMcModelResult = asianMcModel.simulate(simulator,100000, 100);
    cout << "-------------" << endl;
    cout << "Asian Call with 100k paths and 100 steps: " << endl;
    cout << asianMcModelResult;
