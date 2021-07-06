from NumbaLSODA import lsoda_sig, lsoda
import numba as nb
import numpy as np
from scipy.integrate import solve_ivp
import timeit

@nb.cfunc(lsoda_sig)
def f_nb(t, u, du, p):
    k1 = p[0]
    k2 = p[1]
    k3 = p[2]
    du[0] = -k1*u[0]+k3*u[1]*u[2]
    du[1] = k1*u[0]-k2*u[1]**2-k3*u[1]*u[2]
    du[2] = k2*u[1]**2
    
# scipy
@nb.njit
def f_sp(t, u, p):
    k1, k2, k3 = p
    du = np.empty((3,),np.float64)
    du[0] = -k1*u[0]+k3*u[1]*u[2]
    du[1] = k1*u[0]-k2*u[1]**2-k3*u[1]*u[2]
    du[2] = k2*u[1]**2
    return du
    
u0 = np.array([1.0,0.0,0.0])
t_eval = np.linspace(0.0,1.0e5,100)
args_tuple = (np.array([0.04,3e7,1e4]),)
args = np.array([0.04,3e7,1e4])
tspan = [min(t_eval),max(t_eval)]
rtol = 1.0e-8
atol = 1.0e-8

usol_sp = solve_ivp(f_sp, tspan, u0, t_eval=t_eval, args=args_tuple, rtol=rtol, atol=atol, method='LSODA')
funcptr = f_nb.address
usol_nb, success = lsoda(funcptr, u0, t_eval, args, rtol=rtol, atol=atol)

if not np.all(np.isclose(usol_sp.y.T,usol_nb)):
    print("Scipy and NumbaLSODA solutions DO NOT match!!!\n")
    
@nb.njit(boundscheck=False)
def time_nb():
    usol_nb, success = lsoda(funcptr, u0, t_eval, args, rtol=rtol, atol=atol)

def time_sp_LSODA():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='LSODA')

def time_sp_BDF():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='BDF')

def time_sp_Radau():
    usol_sp = solve_ivp(f_sp, tspan, u0, t_eval = t_eval, args=args_tuple, rtol=rtol, atol=atol, method='Radau')
    
# time
time_nb()
time_sp_LSODA()
time_sp_BDF()
time_sp_Radau()
iters = 30
t_nb = timeit.Timer(time_nb).timeit(number=iters)/iters*1e3
t_sp_LSODA = timeit.Timer(time_sp_LSODA).timeit(number=iters)/iters*1e3
t_sp_BDF = timeit.Timer(time_sp_BDF).timeit(number=iters)/iters*1e3
t_sp_Radau = timeit.Timer(time_sp_Radau).timeit(number=iters)/iters*1e3

print("NumbaLSODA time =",'%.3f'%t_nb,'ms') # 0.221 ms
print("Scipy LSODA time =",'%.3f'%t_sp_LSODA,'ms') # 26.483 ms
print("Scipy BDF time =",'%.3f'%t_sp_BDF,'ms') # 158.911 ms
print("Scipy Radau time =",'%.3f'%t_sp_Radau,'ms') # 167.977 ms
