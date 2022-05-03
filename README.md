# numbalsoda

`numbalsoda` is a python wrapper to the LSODA method in [ODEPACK](https://computing.llnl.gov/projects/odepack), which is for solving ordinary differential equation initial value problems. LSODA was originally written in Fortran. `numbalsoda` is a wrapper to a C++ re-write of the original code: https://github.com/dilawar/libsoda 

`numbalsoda` also wraps the `dop853` explicit Runge-Kutta method from this repository: https://github.com/jacobwilliams/dop853

This package is very similar to `scipy.integrate.solve_ivp` ([see here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)), when you set `method = 'LSODA'` or `method = DOP853`. But, `scipy.integrate.solve_ivp` invokes the python interpreter every time step which can be slow. Also, `scipy.integrate.solve_ivp` can not be used within numba jit-compiled python functions. In contrast, `numbalsoda` never invokes the python interpreter during integration and can be used within a numba compiled function which makes `numbalsoda` a lot faster than scipy for most problems, and achieves similar performance to Julia's DifferentialEquations.jl in some cases (see `benchmark` folder).

## Installation
Conda:
```
conda install -c conda-forge numbalsoda
```
Pip:
```
python -m pip install numbalsoda
```

## Basic usage

```python
from numbalsoda import lsoda_sig, lsoda, dop853
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

# integrate with lsoda method
usol, success = lsoda(funcptr, u0, t_eval, data = data)

# integrate with dop853 method
usol1, success1 = dop853(funcptr, u0, t_eval, data = data)

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

## Passing data to the right-hand-side function

In the examples shown above, we passed a an single array of floats to the right-hand-side function:

```python
# ...
data = np.array([1.0])
usol, success = lsoda(funcptr, u0, t_eval, data = data)
```

However, sometimes you might want to pass more data types than just floats. For example, you might want to pass several integers, an array of floats, and an array of integers. One way to achieve this is with generating the `cfunc` using a function like this:

```python
def make_lsoda_func(param1, param2, param3):
    @cfunc(lsoda_sig)
    def rhs(t, x, du, p):
        # Here param1, param2, and param3
        # can be accessed.
        du[0] = param1*t
        # etc...
    return rhs
    
rhs = make_lsoda_func(10.0, 5, 10000)
funcptr = rhs.address
# etc...
```

The only drawback of this approach is if you want to do many successive integrations where the parameters change because it would required re-compiling the `cfunc` between each integration. This could be slow.

But! It is possible to pass arbitrary parameters without re-compiling the `cfunc`, but it is a little tricky. The notebook `passing_data_to_rhs_function.ipynb` gives an example that explains how.
