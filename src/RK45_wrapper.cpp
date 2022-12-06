/**
 * @file RK45_wrapper.cpp
 * @author Thibault Bridel-Bertomeu (re-cae.com)
 * @brief Runge-Kutta-Felhberg method of order 5(4) with adaptive timestepping.
 * @version 1.0
 * @date 2022-12-01
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "RK45.h"

extern "C"
{

/**
 * @brief ODE solver based on Runge-Kutta-Fehlberg 5(4) method.
 * @details This method implements the Runge-Kutta-Fehlberg 5(4) time integrator with adaptive
 * timestepping to solve an ODE.
 * 
 * @param rhs Pointer to the method to compute the right hand side of the ODE.
 * @param neq Number of equations and unknowns
 * @param u0 Initial solution/guess.
 * @param data Extra data to be passed down to the RHS.
 * @param dt0 Initial time step, optional
 * @param t0 Initial time.
 * @param tf Final time
 * @param itf Maximum number of iterations
 * @param usol Vector holding the solution at a given instant.
 * @param tsol Vector holding the instants at which the solution has been evaluated.
 * @param rtol Relative tolerance for the integration.
 * @param atol Absolute tolerance for the integration.
 * @param mxstep Maximum allowed step size.
 * @param success Boolean indicating whether integration was successful or errored.
 */
void rk45_wrapper(
                  void (*rhs)(double t, double *u, double *du, void* data),
                  int neq,
                  double* u0,
                  void* data,
                  double dt0,
                  double t0,
                  double tf,
                  int itf,
                  double* usol,
                  double* tsol,
                  double rtol,
                  double atol,
                  double mxstep,
                  int* success,
                  double* actual_final_time,
                  int* actual_final_iteration
)
{
    std::vector<double> y0(neq, double(0));
    for (auto i = 0; i < neq; i++) {
        y0[i] = u0[i];
    }

    RK45 rk45(rhs, t0, y0, tf, mxstep, rtol, atol, dt0, data);

    for (auto i = 0; i < neq; i++) {
        usol[i + (0) * neq] = rk45.m_y[i];
    }
    tsol[(0)] = rk45.m_t;

    int itnum = 1;
    while (itnum < itf + 1) {
        bool successful_step = rk45.step(data);

        if (not successful_step) {
            *success = 0;
            return;
        }

        for (auto i = 0; i < neq; i++) {
            usol[i + itnum * neq] = rk45.m_y[i];
        }
        tsol[itnum] = rk45.m_t;

        if (rk45.m_direction * (rk45.m_t - rk45.m_t_bound) >= 0) {
            break;
        }

        itnum += 1;
    }

    *actual_final_time = rk45.m_t;
    *actual_final_iteration = itnum;
    *success = 1;

} // end void rk45_wrapper(...)

} // end extern "C"
