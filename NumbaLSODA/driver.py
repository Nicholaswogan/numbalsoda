import ctypes as ct
import numba as nb
from numba import njit, types
import numpy as np
import os
import platform

lsoda_sig = types.void(types.double,
                       types.CPointer(types.double),
                       types.CPointer(types.double),
                       types.CPointer(types.double))

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

if platform.uname()[0] == "Windows":
    name = "liblsoda.dll"
elif platform.uname()[0] == "Linux":
    name = "liblsoda.so"
else:
    name = "liblsoda.dylib"
liblsoda = ct.CDLL(rootdir+name)
lsoda_wrapper = liblsoda.lsoda_wrapper
lsoda_wrapper.argtypes = [ct.c_void_p, ct.c_int, ct.c_void_p, ct.c_void_p,\
                          ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_double,\
                          ct.c_double, ct.c_void_p]
lsoda_wrapper.restype = None

@njit
def lsoda(funcptr, u0, t_eval, data = np.array([0.0], np.float64), rtol = 1.0e-3, atol = 1.0e-6):
    neq = len(u0)
    nt = len(t_eval)
    usol = np.empty((nt,neq),dtype=np.float64)
    success = np.array(999,np.int32)
    lsoda_wrapper(funcptr, neq, u0.ctypes.data, data.ctypes.data, nt,
                  t_eval.ctypes.data, usol.ctypes.data, rtol, atol,
                  success.ctypes.data)
    success_ = True
    if success != 1:
        success_ = False
    return usol, success_
    
@nb.extending.intrinsic
def address_as_void_pointer(typingctx, src):
    """ returns a void pointer from a given memory address """
    from numba import types 
    from numba.core import cgutils
    sig = types.voidptr(src)
    def codegen(cgctx, builder, sig, args):
        return builder.inttoptr(args[0], cgutils.voidptr_t)
    return sig, codegen
