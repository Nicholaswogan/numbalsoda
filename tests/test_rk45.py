from numbalsoda import rk45_sig, rk45
from numba import cfunc
import numpy as np

@cfunc(rk45_sig)
def f(t, u, du, p):
    du[0] = u[0]-u[0]*u[1]
    du[1] = u[0]*u[1]-u[1]

funcptr = f.address

def test_rk45():
    u0 = np.array([5.,0.8])
    data = np.array([1.0])
    dt0 = -1.0
    t0 = 0.0
    tf = 10.0
    itf = 1000000
    
    usol, actual_tf, success = rk45(
        funcptr,
        u0,
        dt0,
        t0,
        tf,
        itf,
        rtol=1.e-8,
        atol=1.e-8
    )

    print(actual_tf)
    print(usol.shape)
    print(usol[-20:,:])

    assert success
    assert np.isclose(usol[-1,0], 0.38583246250193476)
    assert np.isclose(usol[-1,1], 4.602012234037773)

if __name__ == "__main__":
    test_rk45()
