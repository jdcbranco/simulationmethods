
#ifndef SIMULATIONMETHODS_SIMULATOR_H
#define SIMULATIONMETHODS_SIMULATOR_H

#include "Path.h"
#include "Random.h"

using namespace std;

/**
 * The job of the Simulator is generate the simulated paths.
 */
class Simulator {
protected:
    Normal &m_Rng;
    bool m_Antithetic;
public:
    Simulator(Normal &rng, bool antithetic = false): m_Rng(rng), m_Antithetic(antithetic) {}
    vector<double> generate(int size) {
        return m_Rng.generate(size);
    }
    bool is_Antithetic() { return m_Antithetic; }

    vector<Path2> simulate(function<double(double,Bump)> mu_model,
                          function<double (double,Bump)> sigma_model,
                          SDESolver sdeSolver,
                          double S0, double T,
                          double bump_size,
                          int simulations,
                          int path_size = 1) {
        vector<Path2> sims;
        for (int i = 0; i < simulations; ++i) {
            if(m_Antithetic) {
                switch(sdeSolver) {
                    case Explicit:
                        sims.push_back(LogNormalPath(mu_model, sigma_model, m_Rng.generate(path_size), S0, T, bump_size, true));
                        sims.push_back(LogNormalPath(mu_model, sigma_model, m_Rng.generate(path_size), S0, T, bump_size, false));
                        continue;
                    case Euler:
                    case Milstein: //We don't actually support Milstein yet
                        sims.push_back(EulerPath(mu_model, sigma_model, m_Rng.generate(path_size), S0, T, bump_size, true));
                        sims.push_back(EulerPath(mu_model, sigma_model, m_Rng.generate(path_size), S0, T, bump_size, false));
                        continue;
                }
            } else {
                sims.push_back(EulerPath(mu_model, sigma_model, m_Rng.generate(path_size), S0, T, bump_size));
            }
        }
        return sims;
    }

    vector<Path> simulate(ModelParams &model, int simulations, int path_size = 1) {
        vector<Path> sims;
        for (int i = 0; i < simulations; ++i) {
            if(m_Antithetic) {
                sims.push_back(Path(model, m_Rng.generate(path_size), true));
                sims.push_back(Path(model, m_Rng.generate(path_size), false));
            } else {
                sims.push_back(Path(model, m_Rng.generate(path_size)));
            }
        }
        return sims;
    }
};


#endif //SIMULATIONMETHODS_SIMULATOR_H
