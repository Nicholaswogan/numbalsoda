import ctypes as ct
import numba as nb
from numba import njit, types
import numpy as np
import os
import platform

rk45_sig = types.void(types.double,
                      types.CPointer(types.double),
                      types.CPointer(types.double),
                      types.CPointer(types.double))

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

if platform.uname()[0] == "Windows":
    name = "librk45.dll"
elif platform.uname()[0] == "Linux":
    name = "librk45.so"
else:
    name = "librk45.dylib"
librk45 = ct.CDLL(rootdir+name)
rk45_wrapper = librk45.rk45_wrapper
rk45_wrapper.argtypes = [ct.c_void_p, ct.c_int, ct.c_void_p, ct.c_void_p,\
                         ct.c_double, ct.c_double, ct.c_double, ct.c_int,\
                         ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_void_p,\
                         ct.c_void_p, ct.c_void_p]
rk45_wrapper.restype = None

@njit
def rk45(funcptr, u0, dt0, t0, tf, itf, data = np.array([0.0], np.float64), rtol = 1.0e-3, atol = 1.0e-6, mxstep = 10000.0):
    neq = len(u0)
    usol = np.empty((itf+1, neq),dtype=np.float64)
    success = np.array([999], np.int32)
    actual_final_time = np.array([-1.0], dtype=np.float64)
    actual_final_iteration = np.array([-1], dtype=np.int32)
    rk45_wrapper(funcptr, neq, u0.ctypes.data, data.ctypes.data, dt0, t0, tf, itf,
                 usol.ctypes.data, rtol, atol, mxstep, success.ctypes.data,
                 actual_final_time.ctypes.data, actual_final_iteration.ctypes.data)
    success_ = True
    if success != 1:
        success_ = False
    return usol[:actual_final_iteration[0],:], actual_final_time[0], success_
    
@nb.extending.intrinsic
def address_as_void_pointer(typingctx, src):
    """ returns a void pointer from a given memory address """
    from numba import types 
    from numba.core import cgutils
    sig = types.voidptr(src)
    def codegen(cgctx, builder, sig, args):
        return builder.inttoptr(args[0], cgutils.voidptr_t)
    return sig, codegen
