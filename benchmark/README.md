# Benchmarks

This is a comparison `NumbaLSODA`, `scipy.integrate.solve_ivp` and Julia's DifferentialEquations.jl. I look at one non-stiff problem (lorenz) and one stiff problem (rober). Below are the results.

For larger problems (more ODEs), there will be smaller differences between NumbaLSODA and Scipy because Scipy will spend relatively less time doing operations in python (which is slow), and will spend more time in optimized machine code (which is fast).

Also, Scipy might be faster than NumbaLSODA for stiff problems with > ~500 ODEs because Scipy uses highly optimized BLAS and LAPACK routines for matrix inversion (intel MKL or OpenBLAS), while NumbaLSODA uses home-brewed and less optimized BLAS routines.

## Non-stiff (Lorenz)

|              | Time        | Relative to NumbaLSODA |
|--------------|-------------|------------------------|
| NumbaLSODA   | 3.421 ms    | 1x                     |
| Scipy LSODA  | 595.029 ms  | 174x                   |
| Scipy RK45   | 1744.073 ms | 510x                   |
| Scipy DOP853 | 1052.442 ms | 308x                   |
| Julia Tsit5  | 2.612 ms    | 0.76x                  |
| Julia Vern8  | 1.803 ms    | 0.53x                  |

## Stiff (Rober)

|                    | Time       | Relative to NumbaLSODA |
|--------------------|------------|------------------------|
| NumbaLSODA         | 0.221 ms   | 1x                     |
| Scipy LSODA        | 26.483 ms  | 120x                   |
| Scipy BDF          | 158.911 ms | 719x                   |
| Scipy Radau        | 167.977 ms | 760x                   |
| Julia TRBDF2       | 5.765 ms   | 26.08x                 |
| Julia CVODE_BDF    | 1.115 ms   | 5.045x                 |
| Julia LSODA        | 0.264 ms   | 1.194x                 |
| Julia Rodas5       | 0.762 ms   | 3.45x                  |
| Julia Rodas5 + StaticArrays | 0.173 ms   | 0.78x                 |
| Julia Rodas5 + StaticArrays + Analytical Jacobian       | 0.113 ms   | 0.511x                 |
