/**
 * @file RK_merson_wrapper.cpp
 * @author Thibault Bridel-Bertomeu (re-cae.com)
 * @brief Merson Runge Kutta 5(6) ODE integrator.
 * @version 1.0
 * @date 2022-12-01
 */

#include "RK_merson.h"

#include <array>
#include <vector>

const std::array<double, RK_MERSON_NSTEPS> rk_merson_nodes = {0.0, (1.0/3.0), (1.0/3.0), 0.5, 1.0};
const std::array<double, RK_MERSON_NSTEPS> rk_merson_weights = {(1.0/6.0), 0.0, 0.0, 2.0*(1.0/3.0), (1.0/6.0)};
const std::array<double, RK_MERSON_NSTEPS> rk_merson_starweights = {0.5, 0.0, -3.0*0.5, 2.0, 0.0};
const std::array<std::array<double, RK_MERSON_NSTEPS>, RK_MERSON_NSTEPS> rk_merson_matrix = {
    0.0, 0.0, 0.0, 0.0, 0.0,
    (1.0/3.0), 0.0, 0.0, 0.0, 0.0,
    (1.0/6.0), (1.0/6.0), 0.0, 0.0, 0.0,
    (1.0/8.0), 0.0, 3.0*(1.0/8.0), 0.0, 0.0,
    0.5, 0.0, -3.0*0.5, 2.0, 0.0,
};

extern "C"
{

/**
 * @brief ODE solver based on Merson's Runge-Kutta 5(6) integrator.
 * @details This method implements Merson's Runge-Kutta 5(6) method with adaptive time-stepping
 * to solve an ODE.
 * 
 * @param rhs Pointer to the method to compute the right hand side of the ODE.
 * @param neq Number of equations.
 * @param u0 Initial solution/guess.
 * @param data Extra data to be passed down to the RHS.
 * @param nt Size of the teval vector - see below.
 * @param teval Vector of times whereat to solve the ODE (optional).
 * @param usol Vector holding the solution at a given instant.
 * @param rtol Relative tolerance for the integration.
 * @param atol Absolute tolerance for the integration.
 * @param mxstep Maximum number of time steps taken.
 * @param success Boolean indicating whether integration was successful or errored.
 */
void rk_merson_wrapper(
    void (*rhs)(double t, double *u, double *du, void* data),
    int neq,
    double* u0,
    void* data,
    int nt,
    double* teval,
    double* usol,
    double rtol,
    double atol,
    int mxstep,
    int* success
)
{



}

}
