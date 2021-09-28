
#include "LSODA.h"
#include <iostream>
#include <vector>


extern "C"
{
    
void lsoda_wrapper(void (*rhs)(double t, double *u, double *du, void *data),
                   int neq, double* u0, void* data, int nt, double* teval,
                   double* usol, double rtol, double atol, int* success){
  
  LSODA lsoda;  
  std::vector<double> y;
  std::vector<double> yout;
  int istate = 1;
  y.resize(neq);
  yout.resize(neq);
  *success = 1;
  
  // load in initial conditions
  for (int i = 0; i < neq; i++){
    y[i] = u0[i];
    usol[i] = u0[i];
  }
  double t = teval[0];
  double tout;
  
  for (int i = 1; i < nt; i++){
    if (teval[i] < teval[i-1]){
      *success = 0;
      return;
    }
    
    tout = teval[i];
    lsoda.lsoda_update(rhs, neq, y, yout, &t, tout, &istate, data, rtol, atol);
    
    if (istate <= 0){
      // there is a problem!
      *success = 0;
      return;
    }
    
    // update y for next step
    for (int j = 0; j < neq; j++){
      y[j] = yout[j+1];
    }
    
    // save solution
    for (int j = 0; j < neq; j++){
      usol[j + neq*i] = yout[j+1];
    }
  }
}

}