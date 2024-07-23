// SimulationParams.h
#ifndef SIMULATION_PARAMS_H
#define SIMULATION_PARAMS_H

#include <string>

struct simulation_params {
    int forward_euler_only;
    double celltype;
    double bcl;
    int beats;
    double dtw;
    double time_point;
    double min_dt;
    double max_dt;
    std::string drug_name;
    double* conc;
    int conc_size;
};

// Declare the function to load parameters
simulation_params load_params(const char* filename);

#endif // SIMULATION_PARAMS_H
