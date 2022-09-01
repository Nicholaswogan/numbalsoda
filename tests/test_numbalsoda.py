from numbalsoda import lsoda_sig, lsoda, dop853, solve_ivp
from numba import njit, cfunc
import numpy as np

@cfunc(lsoda_sig)
def f(t, u, du, p):
    du[0] = u[0]-u[0]*u[1]
    du[1] = u[0]*u[1]-u[1]
funcptr = f.address

def test_lsoda():
    u0 = np.array([5.,0.8])
    data = np.array([1.0])
    t_eval = np.linspace(0.0,10.0,11)

    usol, success = lsoda(funcptr, u0, t_eval,rtol=1.e-8,atol=1.e-8)

    assert success
    assert np.isclose(usol[-1,0], 0.38583246250193476)
    assert np.isclose(usol[-1,1], 4.602012234037773)

def test_dop853():
    u0 = np.array([5.,0.8])
    data = np.array([1.0])
    t_eval = np.linspace(0.0,10.0,11)

    usol, success = dop853(funcptr, u0, t_eval,rtol=1.e-8,atol=1.e-8)

    assert success
    assert np.isclose(usol[-1,0], 0.38583246250193476)
    assert np.isclose(usol[-1,1], 4.602012234037773)

def test_solve_ivp_1():
    u0 = np.array([5.,0.8])
    data = np.array([1.0])
    t_eval = np.linspace(0.0,10.0,11)
    t_span = np.array([min(t_eval),max(t_eval)])

    sol = solve_ivp(funcptr, t_span, u0, t_eval = t_eval, rtol=1.e-8, atol=1.e-8)

    assert sol.success
    assert np.isclose(sol.y[-1,0], 0.38583246250193476)
    assert np.isclose(sol.y[-1,1], 4.602012234037773)

def test_solve_ivp_2():
    u0 = np.array([5.,0.8])
    data = np.array([1.0])
    t_eval = np.linspace(0.0,10.0,11)
    t_span = np.array([min(t_eval),max(t_eval)])

    sol = solve_ivp(funcptr, t_span, u0, rtol=1.e-8, atol=1.e-8)

    assert sol.success
    assert np.isclose(sol.y[-1,0], 0.38583246250193476)
    assert np.isclose(sol.y[-1,1], 4.602012234037773)


if __name__ == "__main__":
    test_lsoda()
    test_dop853()
    test_solve_ivp_1()
    test_solve_ivp_2()