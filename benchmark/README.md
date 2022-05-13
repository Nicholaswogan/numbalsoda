# Benchmarks

This is a comparison `numbalsoda`, `scipy.integrate.solve_ivp` and Julia's DifferentialEquations.jl. I look at one non-stiff problem (lorenz) and one stiff problem (rober). Below are the results.

For larger problems (more ODEs), there will be smaller differences between numbalsoda and Scipy because Scipy will spend relatively less time doing operations in python (which is slow), and will spend more time in optimized machine code (which is fast).

Also, Scipy might be faster than numbalsoda for stiff problems with > ~500 ODEs because Scipy uses highly optimized BLAS and LAPACK routines for matrix inversion (intel MKL or OpenBLAS), while numbalsoda uses home-brewed and less optimized BLAS routines.

## Non-stiff (Lorenz)

|                            | Time       | Relative to numbalsoda |
| -------------------------- | ---------- | ---------------------- |
| numbalsoda lsoda           | 2.594 ms   | 1x                     |
| numbalsoda dop853          | 0.700 ms   | 0.27x                  |
| Scipy LSODA                | 109.349 ms | 42x                    |
| Scipy RK45                 | 292.826 ms | 113x                   |
| Scipy DOP853               | 180.997 ms | 70x                    |
| Julia Tsit5 + StaticArrays | 1.597 ms   | 0.62x                  |
| Julia Vern8 + StaticArrays | 0.899 ms   | 0.35x                  |

## Stiff (Rober)

|                             | Time     | Relative to numbalsoda |
| --------------------------- | -------- | ---------------------- |
| numbalsoda                  | 0.128 ms | 1x                     |
| Scipy LSODA                 | 4.805 ms | 38x                    |
| Scipy BDF                   | 27.25 ms | 213x                   |
| Scipy Radau                 | 28.34 ms | 221x                   |
| Julia TRBDF2                | 1.524 ms | 11.9x                  |
| Julia CVODE_BDF             | 0.466 ms | 3.6x                   |
| Julia Rodas5                | 0.155 ms | 1.2x                   |
| Julia Rodas5 + StaticArrays | 0.053 ms | 0.41x                  |
