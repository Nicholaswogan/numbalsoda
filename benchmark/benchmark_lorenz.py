import numpy as np
from numbalsoda import lsoda_sig, lsoda, dop853
from scipy.integrate import solve_ivp
import timeit
import numba as nb

# numbalsoda
@nb.cfunc(lsoda_sig,boundscheck=False)
def f_nb(t, u_, du_, p_):
    u = nb.carray(u_, (3,))
    p = nb.carray(p_, (3,))
    sigma, rho, beta = p
    x, y, z = u
    du_[0] = sigma * (y - x)
    du_[1] = x * (rho - z) - y
    du_[2] = x * y - beta * z
funcptr = f_nb.address

# scipy
@nb.njit
def f_sp(t, u, sigma, rho, beta):
    x, y, z = u
    return np.array([sigma * (y - x), x * (rho - z) - y, x * y - beta * z])
    
u0 = np.array([1.0,0.0,0.0])
tspan = np.array([0., 100.])
t_eval = np.linspace(tspan[0], tspan[1], 1001)
args = np.array([10.0,28.0,8./3.])
args_tuple = tuple(args)
rtol = 1.0e-8
atol = 1.0e-8

usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='LSODA')
usol_nb, success = lsoda(funcptr, u0, t_eval, args, rtol=rtol, atol=atol) 

@nb.njit(boundscheck=False)
def time_nb():
    usol_nb, success = lsoda(funcptr, u0, t_eval, args, rtol=rtol, atol=atol)  
    
@nb.njit(boundscheck=False)
def time_dop853():
    usol_nb, success = dop853(funcptr, u0, t_eval, args, rtol=rtol, atol=atol)  

def time_sp_LSODA():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='LSODA')

def time_sp_RK45():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='RK45')

def time_sp_DOP853():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='DOP853')

# time
time_nb()
time_dop853()
time_sp_LSODA()
time_sp_RK45()
time_sp_DOP853()
iters = 10
t_nb = timeit.Timer(time_nb).timeit(number=iters)/iters*1e3
t_dop853 = timeit.Timer(time_dop853).timeit(number=iters)/iters*1e3
t_sp_LSODA = timeit.Timer(time_sp_LSODA).timeit(number=iters)/iters*1e3
t_sp_RK45 = timeit.Timer(time_sp_RK45).timeit(number=iters)/iters*1e3
t_sp_DOP853 = timeit.Timer(time_sp_DOP853).timeit(number=iters)/iters*1e3
print("numbalsoda lsoda time =",'%.3f'%t_nb,'ms')
print("numbalsoda dop853 time =",'%.3f'%t_dop853,'ms')
print("Scipy LSODA time =",'%.3f'%t_sp_LSODA,'ms')
print("Scipy RK45 time =",'%.3f'%t_sp_RK45,'ms')
print("Scipy DOP853 time =",'%.3f'%t_sp_DOP853,'ms')

