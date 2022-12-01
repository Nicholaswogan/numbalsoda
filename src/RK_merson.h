/**
 * @file RK_merson.h
 * @author Thibault Bridel-Bertomeu (re-cae.com)
 * @brief Merson Runge Kutta 5(6) ODE integrator utilities.
 * @version 1.0
 * @date 2022-12-01
 */

#ifndef RK_MERSON_H
#define RK_MERSON_H

#include <array>

#define RK_MERSON_NSTEPS 5

extern const std::array<double, RK_MERSON_NSTEPS> rk_merson_nodes;
extern const std::array<double, RK_MERSON_NSTEPS> rk_merson_weights;
extern const std::array<double, RK_MERSON_NSTEPS> rk_merson_starweights;
extern const std::array<std::array<double, RK_MERSON_NSTEPS>, RK_MERSON_NSTEPS> rk_merson_matrix;

#endif
