# NumbaLSODA

`NumbaLSODA` is a python wrapper to the LSODA method in [ODEPACK](https://computing.llnl.gov/projects/odepack), which is for solving ordinary differential equation initial value problems. LSODA was originally written in Fortran. `NumbaLSODA` is a wrapper to a C++ re-write of the original code: https://github.com/dilawar/libsoda 

This package is very similar to `scipy.integrate.solve_ivp` ([see here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)), when you set `method = 'LSODA'`. But, `scipy.integrate.solve_ivp` invokes the python interpreter every time step which can be slow. Also, `scipy.integrate.solve_ivp` can not be used within numba jit-compiled python functions. In contrast, `NumbaLSODA` never invokes the python interpreter during integration and can be used within a numba compiled function which makes `NumbaLSODA` a lot faster than scipy for most problems (see `benchmark` folder).

## Installation
`NumbaLSODA` should work on Windows, Linux, or MacOS. Install with pip:

```
python -m pip install NumbaLSODA
```

## Basic usage

```python
from NumbaLSODA import lsoda_sig, lsoda
from numba import njit, cfunc
import numpy as np

@cfunc(lsoda_sig)
def rhs(t, u, du, p):
    du[0] = u[0]-u[0]*u[1]
    du[1] = u[0]*u[1]-u[1]*p[0]

funcptr = rhs.address # address to ODE function
u0 = np.array([5.,0.8]) # Initial conditions
data = np.array([1.0]) # data you want to pass to rhs (data == p in the rhs).
t_eval = np.linspace(0.0,50.0,1000) # times to evaluate solution

usol, success = lsoda(funcptr, u0, t_eval, data = data)
# usol = solution
# success = True/False
```

The variables `u`, `du` and `p` in the `rhs` function are pointers to an array of floats. Therefore, operations like `np.sum(u)` or `len(u)` will not work. However, you can use the function `nb.carray()` to make a numpy array out of the pointers. For example:

```python
import numba as nb

@cfunc(lsoda_sig)
def rhs(t, u, du, p):
    u_ = nb.carray(u, (2,))
    p_ = nb.carray(p, (1,))
    # ... rest of rhs goes here using u_ and p_
```

Above, `u_` and `p_` are numpy arrays build out of `u` and `p`, and so functions like `np.sum(u_)` will work.

Also, note `lsoda` can be called within a jit-compiled numba function (see below). This makes it much faster than scipy if a program involves many integrations in a row.

```python
@njit
def test():
    usol, success = lsoda(funcptr, u0, t_eval, data = data)
    return usol
usol = test() # this works!

@njit
def test_sp():
    sol = solve_ivp(f_scipy, t_span, u0, t_eval = t_eval, method='LSODA')
    return sol
sol = test_sp() # this does not work :(
```
