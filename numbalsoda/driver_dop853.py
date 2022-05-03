import ctypes as ct
import numba as nb
from numba import njit, types
import numpy as np
import os
import platform

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

if platform.uname()[0] == "Windows":
    name = "libdop853.dll"
elif platform.uname()[0] == "Linux":
    name = "libdop853.so"
else:
    name = "libdop853.dylib"
libdop853 = ct.CDLL(rootdir+name)
dop853_wrapper = libdop853.dop853_wrapper
dop853_wrapper.argtypes = [ct.c_void_p, ct.c_int, ct.c_void_p, ct.c_void_p,\
                          ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_double,\
                          ct.c_double, ct.c_int, ct.c_void_p]
dop853_wrapper.restype = None

@njit
def dop853(funcptr, u0, t_eval, data = np.array([0.0], np.float64), rtol = 1.0e-3, atol = 1.0e-6, mxstep = 10000):
    neq = len(u0)
    nt = len(t_eval)
    usol = np.empty((nt,neq),dtype=np.float64)
    success = np.array(999,np.int32)
    dop853_wrapper(funcptr, neq, u0.ctypes.data, data.ctypes.data, nt,
                  t_eval.ctypes.data, usol.ctypes.data, rtol, atol, mxstep,
                  success.ctypes.data)
    success_ = True
    if success != 1:
        success_ = False
    return usol, success_

