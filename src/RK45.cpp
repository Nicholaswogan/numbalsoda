/**
 * @file RK45.cpp
 * @author Thibault Bridel-Bertomeu (re-cae.com)
 * @brief Runge-Kutta-Felhberg method of order 5(4); class definition;
 * @version 1.0
 * @date 2022-12-01
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "RK45.h"

/**
 * @brief Construct a new RK45::RK45 object
 * 
 * @param fun Callable, right-hand side of the system.
 * @param t0 Initial time.
 * @param y0 Initial state.
 * @param t_bound Boundary time, the integration won't continue beyond it.
 * @param max_step Maximum allowed step size; default is infinity, i.e. the step size
 * is not bounded and determined solely by the solver.
 * @param rtol Relative tolerance.
 * @param atol Absolute tolerance.
 * @param first_step Initial step size; default is "None" which means that the algorithm
 * should choose.
 */
RK45::RK45(RK45_ODE_SYSTEM_TYPE rhs_handle, double t0, std::vector<double> &y0, double tf,
           int largest_step,
           double _rtol,
           double _atol,
           std::optional<double> first_step)
{
    _m_status = "running";
    m_fun = rhs_handle;

    _m_t_old = std::nullopt;
    m_t = t0;
    m_t_bound = tf;
    m_direction = (m_t_bound != t0) ? sgn(m_t_bound - t0) : 1;

    if (_rtol < 100 * std::numeric_limits<double>::epsilon()) {
        std::cout << "Warning: '_rtol' is too small. Setting to max(rtol, 100*EPS)." << std::endl;
        m_rtol = std::max(_rtol, 100 * std::numeric_limits<double>::epsilon());
    }
    if (_atol < 0) {
        throw std::invalid_argument("'_atol' must be positive.");
    }
    m_atol = _atol;

    m_y = y0;
    _m_y_old = std::nullopt;
    m_n_eqns = m_y.size();
    m_fun(m_t, m_y.data(), m_f.data(), NULL);


    if (largest_step <= int(0))
        throw std::invalid_argument("'largest_step' must be positive.");
    m_max_step = largest_step;
    if (first_step.has_value()) {
        if (first_step.value() <= 0)
            throw std::invalid_argument("'first_step' must be positive.");
        if (first_step.value() > std::abs(tf - t0))
            throw std::invalid_argument("'first_step' exceeds bounds.");
        m_h_abs = first_step.value();
    } else {
        m_h_abs = select_initial_step(m_fun, m_t, m_y, m_f, m_direction, _m_error_estimation_order, m_rtol, m_atol);
    }
    _m_h_previous = std::nullopt;

    _m_K.reserve(_m_n_stages + 1);
    for (auto i = 0; i < _m_n_stages + 1; i++)
        _m_K[i].reserve(m_n_eqns);
}

/**
 * @brief Estimate the current integration error.
 * 
 * @param h Current step size.
 * @return std::vector<double> Current integration error per equation.
 */
std::vector<double> RK45::_estimate_error(double h)
{
    std::vector<double> dot_K_E(m_n_eqns);
    dgemv("T", _m_n_stages + 1, m_n_eqns, 1.0, _m_K, _m_n_stages + 1, _m_E, 1, 0.0, dot_K_E, 1);

    for (auto i = 0; i < dot_K_E.size(); i++) {
        dot_K_E[i] *= h;
    }

    return dot_K_E;
}

/**
 * @brief Norm of the current integration error.
 * 
 * @param h Current step size.
 * @param scale Scaling factor for the current integration error per equation.
 * @return double Norm of the current integration error, each equation rescaled by the
 * ad hoc factor based on user specified tolerances.
 */
double RK45::_estimate_error_norm(double h, std::vector<double> scale)
{
    std::vector<double> itg_error = _estimate_error(h);
    double norm = double(0);
    for (auto i = 0; i < scale.size(); i++) {
        norm += (itg_error[i] / scale[i]) * (itg_error[i] / scale[i]);
    }

    return std::sqrt(norm);
}

/**
 * @brief Computations for one integration step.
 *
 * @return bool Whether the step was successful.
 */
bool RK45::step()
{
    double min_step = double(10) * std::numeric_limits<double>::epsilon();
    double h_abs;
    if (m_h_abs > m_max_step) {
        h_abs = m_max_step;
    } else if (m_h_abs < min_step) {
        h_abs = min_step;
    } else {
        h_abs = m_h_abs;
    }

    bool step_accepted = false;
    bool step_rejected = false;

    while (not step_accepted) {
        if (h_abs < min_step) {
            return false;
        }

        double h = h_abs * m_direction;
        double t_new = m_t + dt;
        if (m_direction * (t_new - m_t_bound) > double(0)) {
            t_new = m_t_bound;
        }
        h = t_new - m_t;
        h_abs = std::abs(h);
    }
}
