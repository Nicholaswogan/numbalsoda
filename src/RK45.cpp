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
 * @brief Empirically select a good initial step.
 * @details The algorithm is described in [1].
 * 
 * References:
 * - E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary Differential
 * Equations I: Nonstiff problems", Sec. II.4.
 * 
 * @param fun Callable, right-hand side of the system.
 * @param t0 Initial value of the independant variable.
 * @param y0 Initial value of the dependant variable.
 * @param f0 Initial value of the derivatives.
 * @param direction Integration direction.
 * @param order Error estimation order, i.e. the error controlled by the
 * algorithm is proportional to 'step_size ** (order + 1)'.
 * @param rtol Desired relative tolerance.
 * @param atol Desired absolute tolerance.
 * @return double Absolute value of the suggested initial step.
 */
double select_initial_step(
    RK45_ODE_SYSTEM_TYPE fun,
    double t0,
    std::vector<double> &y0,
    std::vector<double> &f0,
    int direction,
    int order,
    double rtol,
    double atol
)
{
    int n_eqns = y0.size();
    if (n_eqns == 0)
        return std::numeric_limits<double>::max();

    std::vector<double> scale(n_eqns);
    double d0 = 0.0, d1 = 0.0;
    for (auto i = 0; i < n_eqns; i++) {
        scale[i] = atol + std::abs(y0[i]) * rtol;
        d0 += (y0[i] / scale[i]) * (y0[i] / scale[i]);
        d1 += (f0[i] / scale[i]) * (f0[i] / scale[i]);
    }
    d0 = std::sqrt(d0);
    d1 = std::sqrt(d1);

    double h0;
    if ((d0 < 1e-5) || (d1 < 1e-5)) {
        h0 = 1e-6;
    } else {
        h0 = 0.01 * d0 / d1;
    }

    std::vector<double> y1(n_eqns), f1(n_eqns);
    for (auto i = 0; i < n_eqns; i++) {
        y1[i] = y0[i] + h0 * direction * f0[i];
    }
    fun(t0 + h0 * direction, y1.data(), f1.data(), NULL);
    double d2 = 0.0;
    for (auto i = 0; i < n_eqns; i++) {
        d2 += (f1[i] - f0[i]) / scale[i];
    }
    d2 = std::sqrt(d2) / h0;

    double h1;
    if ((d1 <= 1e-15) && (d2 <= 1e-15)) {
        h1 = std::max(1e-06, h0 * 1e-03);
    } else {
        h1 = std::pow((0.01 / std::max(d1, d2)), 1.0 / (order + 1));
    }

   return std::min(100 * h0, h1);
}

/**
 * @brief Perform a single Runge-Kutta step.
 * @details This method is integrator-agnostic and is valid for all explicit Runge-Kutta
 * built with a Butcher tableau.
 * 
 * @param fun Callable, right-hand side of the system.
 * @param t Current value of the independent variable.
 * @param y Current value of the dependent variable.
 * @param f Current value of the derivative.
 * @param h Step to use.
 * @param A Coefficients for combining previous RK stages to compute the next stage.
 * @param B Coefficients for combining RK stages for computing the final prediction.
 * @param C Coefficients for incrementing time for consecutive RK stages.
 * @param K Storage array for putting RK stages here. Stages are stored in rows. The last
 * row is a linear combination of the previous rows with coefficients.
 * @return std::array<std::vector<double>, 2> Solution at t + h and derivative thereof.
 */
std::array<std::vector<double>, 2> rk_step(
    RK45_ODE_SYSTEM_TYPE fun,
    double const t,
    std::vector<double> const& y,
    std::vector<double> const& f,
    double const h,
    std::vector<std::vector<double>> const& A,
    std::vector<double> const& B,
    std::vector<double> const& C,
    std::vector<std::vector<double>>& K
)
{
    int const n_eqns = y.size();
    int const n_stages = K.size() - 1;
    
    for (auto i = 0; i < n_eqns; i++) {
        K[0][i] = f[i];
    }
    for (auto s = 1; s < n_stages; s++) {
        std::vector<double> dy(n_eqns, double(0));
        for (auto iss = 0; iss < s; iss++) {
            for (auto i = 0; i < n_eqns; i++) {
                dy[i] += K[iss][i] * A[s][iss] * h;
            }
        }
        for (auto i = 0; i < n_eqns; i++) {
            dy[i] += y[i];
        }
        fun(t + C[s] * h, dy.data(), K[s].data(), NULL);
    }

    std::vector<double> y_new(n_eqns, double(0)), f_new(n_eqns);
    for (auto s = 0; s < n_stages; s++) {
        for (auto i = 0; i < n_eqns; i++) {
            y_new[i] += K[s][i] * B[s] * h;
        }
    }
    for (auto i = 0; i < n_eqns; i++) {
        y_new[i] += y[i];
    }
    fun(t + h, y_new.data(), f_new.data(), NULL);
    for (auto i = 0; i < n_eqns; i++) {
        K[n_stages][i] = f_new[i];
    }

    std::array<std::vector<double>, 2> result{y_new, f_new};
    return result;
}

/**
 * @brief Performs y = alpha*A*x + beta*y or y = alpha*A^T*x + beta*y
 * 
 * @param trans Specifies the operation to be performed. If equals "N" then
 * y = alpha * A * x + beta * y, if equals "T" then y = alpha * A^T * x + beta * y.
 * @param m Number of rows of the matrix A.
 * @param n Number of columns of the matrix A.
 * @param alpha Scalar alpha.
 * @param a m x n matrix of coefficients A.
 * @param lda First dimension of A.
 * @param x Vector X.
 * @param incx Increment for the elements of X.
 * @param beta Scalar beta.
 * @param y Vector Y; will be overwritten by the updated vector Y.
 * @param incy Increment for the elements of Y.
 */
void dgemv(const std::string trans, const size_t m, const size_t n,
           const double alpha, const std::vector<std::vector<double>> &a,
           const int lda, const std::vector<double> &x, const int incx,
           const double beta, std::vector<double> &y, const int incy)
{

    int lenx, leny;
    if (trans == "N") {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    int kx, ky;
    if (incx > 0) {
        kx = 0;
    } else {
        kx = - (lenx - 1) * incx;
    }
    if (incy > 0) {
        ky = 0;
    } else {
        ky = - (leny - 1) * incy;
    }

    if (std::abs(beta - double(1)) > 0) {
        if (incy == 1) {
            if (std::abs(beta) <= double(0)) {
                std::fill(y.begin(), y.end(), double(0));
            } else {
                for (auto i = 0; i < leny; i++)
                    y[i] *= beta;
            }
        } else {
            int iy = ky;
            if (std::abs(beta) <= double(0)) {
                for (auto i = 0; i < leny; i++) {
                    y[iy] = double(0);
                    iy += incy;
                }
            } else {
                for (auto i = 0; i < leny; i++) {
                    y[iy] *= beta;
                    iy += incy;
                }
            }
        }
    }

    // Form y = alpha * A * x + beta * y
    if (trans == "N") {
        int jx = kx;
        if (incy == 1) {
            for (auto j = 0; j < n; j++) {
                double temp = alpha * x[jx];
                for (auto i = 0; i < m; i++) {
                    y[i] = y[i] + temp * a[i][j];
                }
                jx += incx;
            }
        } else {
            for (auto j = 0; j < n; j++) {
                double temp = alpha * x[jx];
                int iy = ky;
                for (auto i = 0; i < m; i++) {
                    y[iy] = y[iy] + temp * a[i][j];
                    iy += incy;
                }
                jx += incx;
            }
        }
    // Form y = alpha * A^T * x + beta * y
    } else {
        int jy = ky;
        if (incx == 1) {
            for (auto j = 1; j < n; j++) {
                double temp = double(0);
                for(auto i = 0; i < m; i++) {
                    temp += a[i][j] * x[i];
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            for (auto j = 1; j < n; j++) {
                double temp = double(0);
                int ix = kx;
                for (auto i = 0; i < m; i++) {
                    temp += a[i][j] * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        }
    }
}

/**
 * @brief Construct a new RK45::RK45 object.
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
           double largest_step,
           double _rtol,
           double _atol,
           double first_step)
{
    _m_status = "running";
    m_fun = rhs_handle;

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
    m_n_eqns = m_y.size();
    m_f.resize(m_n_eqns);
    m_fun(m_t, m_y.data(), m_f.data(), NULL);

    if (largest_step <= double(0))
        throw std::invalid_argument("'largest_step' must be positive.");
    m_max_step = largest_step;
    if (first_step > 0) {
        if (first_step > std::abs(tf - t0))
            throw std::invalid_argument("'first_step' exceeds bounds.");
        m_h_abs = first_step;
    } else {
        m_h_abs = select_initial_step(m_fun, m_t, m_y, m_f, m_direction, _m_error_estimation_order, m_rtol, m_atol);
    }
    _m_h_previous = -1.0;

    _m_K.resize(_m_n_stages + 1);
    for (auto i = 0; i < _m_n_stages + 1; i++)
        _m_K[i].resize(m_n_eqns);
}

/**
 * @brief Destroy the RK45::RK45 object.
 */
RK45::~RK45() {};

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

    double h, t_new;
    std::vector<double> y_new(m_n_eqns,  double(0)), f_new(m_n_eqns,  double(0));
    while (not step_accepted) {
        if (h_abs < min_step) {
            return false;
        }

        h = h_abs * m_direction;
        t_new = m_t + h;
        if (m_direction * (t_new - m_t_bound) > double(0)) {
            t_new = m_t_bound;
        }
        h = t_new - m_t;
        h_abs = std::abs(h);

        auto y_and_f_new = rk_step(m_fun, m_t, m_y, m_f, h, _m_A, _m_B, _m_C, _m_K);
        y_new = y_and_f_new[0];
        f_new = y_and_f_new[1];

        std::vector<double> scale(m_n_eqns);
        for (auto i = 0; i < m_n_eqns; i++) {
            scale[i] = m_atol + std::max(std::abs(m_y[i]), std::abs(y_new[i])) * m_rtol;
        }
        double error_norm = _estimate_error_norm(h, scale);

        double factor = double(0);
        if (error_norm < 1) {
            if (std::abs(error_norm) < std::numeric_limits<double>::epsilon()) {
                factor = MAX_FACTOR;
            } else {
                factor = std::min(MAX_FACTOR, SAFETY * std::pow(error_norm, _m_error_exponent));
            }

            if (step_rejected) {
                factor = std::min(double(1), factor);
            }

            h_abs *= factor;

            step_accepted = true;
        } else {
            h_abs *= std::max(MIN_FACTOR, SAFETY * std::pow(error_norm, _m_error_exponent));
            step_rejected = true;
        }
    }

    _m_h_previous = h;
    m_t = t_new;
    for (auto i = 0; i < m_n_eqns; i++) {
        m_y[i] = y_new[i];
    }
    m_h_abs = h_abs;
    for (auto i = 0; i < m_n_eqns; i++) {
        m_f[i] = f_new[i];
    }

    return true;
}
