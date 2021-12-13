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
| Julia Tsit5  | 1.340 ms    | 0.39x                  |
| Julia Vern8  | 1.803 ms    | 0.21x                  |

## Stiff (Rober)

|                    | Time       | Relative to NumbaLSODA |
|--------------------|------------|------------------------|
| NumbaLSODA         | 0.221 ms   | 1x                     |
| Scipy LSODA        | 26.483 ms  | 120x                   |
| Scipy BDF          | 158.911 ms | 719x                   |
| Scipy Radau        | 167.977 ms | 760x                   |
| Julia TRBDF2       | 2.202 ms   | 9.96x                 |
| Julia CVODE_BDF    | 0.628 ms   | 2.84x                 |
| Julia LSODA        | 0.254 ms   | 1.15x                 |
| Julia Rodas5       | 0.227 ms   | 1.03x                  |
| Julia Rodas5 + StaticArrays | 0.0423 ms   | 0.19x                 |
| Julia Rodas5 + StaticArrays + Analytical Jacobian       | 0.0420 ms ms   | 0.19x                 |
