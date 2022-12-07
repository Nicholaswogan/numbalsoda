/**
 * @file RK45.h
 * @author Thibault Bridel-Bertomeu (re-cae.com)
 * @brief Runge-Kutta-Felhberg method of order 5(4); class declaration.
 * @version 1.0
 * @date 2022-12-01
 */

#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Multiply steps computed from asymptotic behavior of errors by this.
 */
#define SAFETY 0.9

/**
 * @brief Minimum allowed decrease in a step size.
 */
#define MIN_FACTOR 0.2

/**
 * @brief Maximum allowed increased in a step size.
 */
#define MAX_FACTOR 10.0

/**
 * @brief Type definition of RK45 ode system.
 */
typedef void (*RK45_ODE_SYSTEM_TYPE)(double t, double *y, double *dydt, void*);

/**
 * @brief Runge-Kutta-Fehlberg 5(4) explicit integration method.
 * @details This uses the Dormand-Prince pair of formulas [1]. The error is controlled
 * assuming accuracy of the fourth-order method accuracy, but steps are taken
 * using the fifth-order accurate formula (local extrapolation is done).
 * A quartic interpolation polynomial is used for the dense output [2].
 * 
 * References:
 * 1. J. R. Dormand, P. J. Prince, "A family of embedded Runge-Kutta
 *   formulae", Journal of Computational and Applied Mathematics, Vol. 6,
 *   No. 1, pp. 19-26, 1980.
 * 2. L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
 *   of Computation,, Vol. 46, No. 173, pp. 135-150, 1986.
 */
class RK45 {

public:
    RK45(RK45_ODE_SYSTEM_TYPE fun, double t0, std::vector<double> &y0, double tf,
         double largest_step=std::numeric_limits<double>::max(),
         double _rtol=1e-3,
         double _atol=1e-6,
         double first_step=double(-1),
         void *const rhs_args=NULL);
    ~RK45();

    RK45_ODE_SYSTEM_TYPE m_fun;
    int m_n_eqns;
    double m_max_step;
    double m_t_bound;
    double m_rtol;
    double m_atol;
    
    std::vector<double> m_y;
    std::vector<double> m_f;
    double m_h_abs;
    double m_t;
    int m_direction;

    bool step(void *const rhs_args);

private:
    static int const _m_order = 5;
    static int const _m_error_estimation_order = 4;
    static constexpr double const _m_error_exponent = -1.0 / (_m_error_estimation_order + 1.0);
    static int const _m_n_stages = 6;

    std::vector<double> const _m_C {0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0}; // size _n_stages
    std::vector<std::vector<double>> const _m_A {
        { 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 1.0/5.0, 0.0, 0.0, 0.0, 0.0 },
        { 3.0/40.0, 9.0/40.0, 0.0, 0.0, 0.0 },
        { 44.0/45.0, -56.0/15.0, 32.0/9.0, 0.0, 0.0 },
        { 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0 },
        { 9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0 },
    }; // size _n_stages x _n_stages
    std::vector<double> const _m_B {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}; // size _n_stages
    std::vector<double> const _m_E {-71.0/57600.0, 0.0, 71.0/16695.0, -71.0/1920.0, 17253.0/339200.0, -22.0/525.0, 1.0/40.0}; // size _n_stages + 1

    double _m_h_previous;
    std::vector<std::vector<double>> _m_K;
    std::string _m_status;

    std::vector<double> _estimate_error(double h);
    double _estimate_error_norm(double h, std::vector<double> scale);
};


/**
 * @brief Implements signum (-1, 0, 1).
 * 
 * @tparam T Type of the number whose sign to evaluate.
 * @param val Number whose sign to evaluate.
 * @return int -1 if the number is negative, 1 if the number is positive, 0 elsehow.
 */
template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
