# Benchmarks

This is a comparison NumbaLSODA, scipy.integrate.solveivp and Julia's DifferentialEquations.jl. I look at one non-stiff problem (lorenz) and one stiff problem (rober). Below are the results.

For larger problems (more ODEs), there will be much smaller differences between NumbaLSODA and Scipy. Scipy will be faster for larger problems because relatively less time is spend setting up the problem in python (slow), and more time is spent within algorithms which are in some cases machine code (fast).

## Non-stiff (Lorenz)

|              | Time        | Relative to NumbaLSODA |
|--------------|-------------|------------------------|
| NumbaLSODA   | 3.421 ms    | 1x                     |
| Scipy LSODA  | 595.029 ms  | 174x                   |
| Scipy RK45   | 1744.073 ms | 510x                   |
| Scipy DOP853 | 1052.442 ms | 308x                   |
| Julia Tsit5  | 2.612 ms    | 0.76x                  |

## Stiff (Rober)

|                    | Time       | Relative to NumbaLSODA |
|--------------------|------------|------------------------|
| NumbaLSODA         | 0.221 ms   | 1x                     |
| Scipy LSODA        | 26.483 ms  | 120x                   |
| Scipy BDF          | 158.911 ms | 719x                   |
| Scipy Radau        | 167.977 ms | 760x                   |
| Julia Rosenbrock32 | 1.035 ms   | 4.68x                  |
| Julia TRBDF2       | 5.765 ms   | 26.08x                 |
| Julia CVODE_BDF    | 1.115 ms   | 5.045x                 |
| Julia LSODA        | 0.264 ms   | 1.194x                 |
