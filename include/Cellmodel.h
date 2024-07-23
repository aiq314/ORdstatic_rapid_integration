#ifndef CELLMODEL_H
#define CELLMODEL_H


class Cellmodel
{
protected:
  Cellmodel(){}
public:
  bool isEctopic;
  bool isS1;
  unsigned int algebraic_size;
  unsigned int constants_size;
  unsigned int states_size;
  double *ALGEBRAIC;
  double *CONSTANTS;
  double *RATES;
  double *STATES;
  virtual ~Cellmodel() {}
  virtual void initConsts() = 0;
  virtual void initConsts(double type){}
  virtual void initConsts (double type,
                           bool is_dutta){}
  virtual void initConsts (double type,
                           double conc,
                           double *hill){}
  virtual void initConsts (double type,
                           double conc,
                           const double *hill,
                           const double *herg){}
  virtual void initConsts (double type,
                           double conc,
                           double *hill,
                           bool is_dutta){}
  virtual void computeRates(double TIME,
                            double *CONSTANTS,
                            double *RATES,
                            double *STATES,
                            double *ALGEBRAIC) = 0;
  virtual void solveAnalytical( double dt) {};
  virtual void solveEuler(double dt){};
  virtual double set_time_step (double TIME,
                                       double time_point,
                                       double min_time_step,
                                       double max_time_step,
                                       double min_dV,
                                       double max_dV,
                                       double* CONSTANTS,
                                       double* RATES,
                                       double* STATES,
                                       double* ALGEBRAIC) {};
  virtual void gaussElimination(double *A, double *b, double *x, int N){};
  virtual void solveRK4(double TIME, double dt) {};
  virtual void mat_vec_multiply(double** a, double* x, double *result, int n) {};
  virtual void solve_rk_hERG(double** a, double* y, double dt, int n ){};
};

#endif // CELLMODEL_H
