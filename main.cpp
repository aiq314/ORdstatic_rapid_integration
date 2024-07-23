#include "Ohara_Rudy_2011.hpp"

#include "omp.h"
#include "simulationparams.h"
#include "drug_data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>

int main(int argc, char* argv[]){
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << "<params.txt> <hill_data.csv>" << std::endl;
    return 1;
  }
  // Load params
  simulation_params params = load_params(argv[1]);
  std::cout << "Running simulation with parameters:" << std::endl;
  std::cout << "forward_euler_only : " << params.forward_euler_only << std::endl;
  std::cout << "celltype : " << params.celltype << std::endl;
  std::cout << "bcl (ms): " << params.bcl << std::endl;
  std::cout << "beats : " << params.beats << std::endl;
  std::cout << "dtw (ms): " << params.dtw << std::endl;
  std::cout << "time_point (ms): " << params.time_point << std::endl;
  std::cout << "min_dt (ms): " << params.min_dt << std::endl;
  std::cout << "max_dt (ms): " << params.max_dt << std::endl;
  std::cout << "drug_name : " << params.drug_name << std::endl;
  std::cout << "conc (nM): ";
  for (int i = 0; i < params.conc_size; i++){
    std::cout<< params.conc[i] << " ";
  }
  std::cout << std::endl;
  // Input for adaptive time step
  double time_point = params.time_point;
  double min_dt = params.min_dt;
  double max_dt = params.max_dt;
  // Load drug data
  hill_data hill = load_hill(argv[2]);
  std::cout << "hill data:" << std::endl;
  for(int i = 0; i < 14; i++){
    std::cout << "Col " << i << ": " << hill.hill[i] << std::endl;
  }
  // Measure time
  clock_t start, end;
  double cpu_time_used;
  for (int conc_id = 0; conc_id < params.conc_size; conc_id++){
    // Input for AP simulation
    double t_max = (double) params.beats * params.bcl;
    double t_curr = 0.0;
    double dt = 0.0;
    double dtw = params.dtw;
    double t_next = t_curr + dtw;
    double next_write_time = t_max - params.bcl;
    double epsilon = 1e-6;
    double conc = params.conc[conc_id];
    // Start of calculations
    int imax = int((t_max - next_write_time) / dtw) + 1;// + ((int(t_max) - int(next_write_time)) % int(dtw) == 0 ? 0 : 1);
    Cellmodel* p_elec;
    std::ofstream vmcheck;
    std:: ostringstream file_name;
    std::string vmcheck_name;
    start = clock();
    p_elec = new Ohara_Rudy_2011();
    p_elec->initConsts(0.0,conc,hill.hill); // drug effects
    p_elec->CONSTANTS[BCL] = params.bcl;
    if (params.forward_euler_only == 1){
      file_name << "vmcheck_";
      file_name << params.drug_name << "_";
      file_name << conc << ".plt";
    } else {
      file_name << "vmcheck_";
      file_name << params.drug_name << "_";
      file_name << conc << "_";
      file_name << params.min_dt << "_";
      file_name << params.max_dt << ".plt";
    }
    vmcheck_name = file_name.str();
    vmcheck.open(vmcheck_name.c_str());
    vmcheck << "Time" << "\t";
    vmcheck << "dt" << "\t";
    vmcheck << "INa" << "\t";
    vmcheck << "INaL" << "\t";
    vmcheck << "Ito" << "\t";
    vmcheck << "ICaL" << "\t";
    vmcheck << "ICaNa" << "\t";
    vmcheck << "ICaK" << "\t";
    vmcheck << "IKr" << "\t";
    vmcheck << "IKs" << "\t";
    vmcheck << "IK1" << "\t";
    vmcheck << "INaCa_i" << "\t";
    vmcheck << "INaCa_ss" << "\t";
    vmcheck << "INaK" << "\t";
    vmcheck << "INab" << "\t";
    vmcheck << "IKb" << "\t";
    vmcheck << "IpCa" << "\t";
    vmcheck << "ICab" << "\t";
    for(int i = 0; i < p_elec->states_size; i++) {
      vmcheck << "STATES[" << i << "]";
      if(i < p_elec->states_size-1) {// Add a tab for all but the last element
        vmcheck << "\t";
      } else {// End the line after the last element
        vmcheck << "\n";
      }
    }
    int iprint = 0;
    while(iprint<imax){
      p_elec->computeRates(t_curr,
                          p_elec->CONSTANTS,
                          p_elec->RATES,
                          p_elec->STATES,
                          p_elec->ALGEBRAIC);
      if (params.forward_euler_only == 1){
        dt = min_dt;
        if (t_curr + dt >= next_write_time) {
          dt = next_write_time - t_curr;
        }
        p_elec->solveEuler(dt);
      } else {
        dt = p_elec->set_time_step(t_curr,
                                  time_point,
                                  0.005,
                                  max_dt,
                                  0.2,
                                  0.8,
                                  p_elec->CONSTANTS,
                                  p_elec->RATES,
                                  p_elec->STATES,
                                  p_elec->ALGEBRAIC);
                if (t_curr + dt >= next_write_time) {
          dt = next_write_time - t_curr;
        }
        p_elec->solveAnalytical(dt);
      }
        
      t_curr = t_curr + dt;
      if (t_curr >= next_write_time){
        vmcheck << iprint * dtw << "\t";
        vmcheck << dt << "\t";
        vmcheck << p_elec->ALGEBRAIC[INa] << "\t";
        vmcheck << p_elec->ALGEBRAIC[INaL] << "\t";
        vmcheck << p_elec->ALGEBRAIC[Ito] << "\t";
        vmcheck << p_elec->ALGEBRAIC[ICaL] << "\t";
        vmcheck << p_elec->ALGEBRAIC[ICaNa] << "\t";
        vmcheck << p_elec->ALGEBRAIC[ICaK] << "\t";
        vmcheck << p_elec->ALGEBRAIC[IKr] << "\t";
        vmcheck << p_elec->ALGEBRAIC[IKs] << "\t";
        vmcheck << p_elec->ALGEBRAIC[IK1] << "\t";
        vmcheck << p_elec->ALGEBRAIC[INaCa_i] << "\t";
        vmcheck << p_elec->ALGEBRAIC[INaCa_ss] << "\t";
        vmcheck << p_elec->ALGEBRAIC[INaK] << "\t";
        vmcheck << p_elec->ALGEBRAIC[INab] << "\t";
        vmcheck << p_elec->ALGEBRAIC[IKb] << "\t";
        vmcheck << p_elec->ALGEBRAIC[IpCa] << "\t";
        vmcheck << p_elec->ALGEBRAIC[ICab] << "\t";
        for(int i = 0; i < p_elec->states_size; i++){
          vmcheck << p_elec->STATES[i];
          if(i < p_elec->states_size-1) {// Add a tab for all but the last element
              vmcheck << "\t";
          } else {// End the line after the last element
          vmcheck << "\n";
          }
        }
        // Increment next_write_time by dtw
        next_write_time += dtw;
        iprint += 1;
      }
    }
    vmcheck.close();
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "Simulation is successfully done\n";
    std::cout << "Computational time: " << cpu_time_used << " seconds\n";
    delete p_elec;
  }
  return 0;
}
