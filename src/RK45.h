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
         int largest_step=std::numeric_limits<int>::max(),
         double _rtol=1e-3,
         double _atol=1e-6,
         std::optional<double> first_step=std::nullopt);
    ~RK45();

    RK45_ODE_SYSTEM_TYPE m_fun;
    int m_n_eqns;
    int m_max_step;
    double m_t_bound;
    double m_rtol;
    double m_atol;
    
    std::vector<double> m_y;
    std::vector<double> m_f;
    double m_h_abs;
    double m_t;
    int m_direction;

    bool step();

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
        { 9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/1865.0 },
    }; // size _n_stages x _n_stages
    std::vector<double> const _m_B {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}; // size _n_stages
    std::vector<double> const _m_E {-71.0/57600.0, 0.0, 71.0/16695.0, -71.0/1920.0, 17253.0/339200.0, -22.0/525.0, 1.0/40.0}; // size _n_stages + 1

    std::optional<double> _m_h_previous;
    std::optional<double> _m_t_old;
    std::optional<std::vector<double>> _m_y_old;
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
