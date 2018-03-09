
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
        return move(sims);
    }
};


#endif //SIMULATIONMETHODS_SIMULATOR_H
