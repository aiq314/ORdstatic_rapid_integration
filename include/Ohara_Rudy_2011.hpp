#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "Cellmodel.h"
#include "enum_Ohara_Rudy_2011.hpp"

#define EULER

class Ohara_Rudy_2011 : public Cellmodel
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void initConsts ();
  void initConsts(double type);
  void initConsts(double type, double conc, double *hill );
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt);
  void solveEuler(double dt);
  void solveRK4(double TIME,double dt);
  double set_time_step(double TIME,
                       double time_point,
                       double min_time_step,
                       double max_time_step,
                       double min_dV,
                       double max_dV,
                       double* CONSTANTS,
                       double* RATES,
                       double* STATES,
                       double* ALGEBRAIC);
private:
  void ___applyDrugEffect(double conc, const double *hill);
  void ___initConsts(double type);
  void ___applyDutta();
};


#endif

