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
#include <vector>

struct RK45 {
    static const int n_stages = 6;
    const std::array<double, n_stages> rk_merson_nodes = {0.0, (1.0/3.0), (1.0/3.0), 0.5, 1.0};
    const std::array<double, n_stages> rk_merson_weights = {(1.0/6.0), 0.0, 0.0, 2.0*(1.0/3.0), (1.0/6.0)};
    const std::array<double, n_stages> rk_merson_starweights = {0.5, 0.0, -3.0*0.5, 2.0, 0.0};
    const std::array<std::array<double, n_stages>, n_stages> rk_merson_matrix = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        (1.0/3.0), 0.0, 0.0, 0.0, 0.0,
        (1.0/6.0), (1.0/6.0), 0.0, 0.0, 0.0,
        (1.0/8.0), 0.0, (3.0/8.0), 0.0, 0.0,
        0.5, 0.0, -1.5, 2.0, 0.0
    };
};

extern "C"
{

/**
 * @brief ODE solver based on Merson's Runge-Kutta 5(6) integrator.
 * @details This method implements Merson's Runge-Kutta 5(6) method with adaptive time-stepping
 * to solve an ODE.
 * 
 * @param rhs Pointer to the method to compute the right hand side of the ODE.
 * @param neq Number of equations and unknowns
 * @param data Extra data to be passed down to the RHS.
 * @param t0 Initial time.
 * @param u0 Initial solution/guess.
 * @param tf Final time
 * @param itf Maximum number of iterations
 * @param rtol Relative tolerance for the integration.
 * @param atol Absolute tolerance for the integration.
 * @param dt0 Initial time step, optional
 * @param usol Vector holding the solution at a given instant.
 * @param success Boolean indicating whether integration was successful or errored.
 */
void rk_merson_wrapper(
                       void (*rhs)(double t, double *u, double *du, void* data),
                       int neq,
                       void* data,
                       double t0,
                       double* u0,
                       double tf,
                       int itf,
                       double rtol,
                       double atol,
                       double dt0,
                       double* usol,
                       int* success
)
{

} // end void rk_merson_wrapper(...)

} // end extern "C"
