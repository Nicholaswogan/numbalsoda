import ctypes as ct
import numba as nb
from numba import types
import numpy as np
import os
import platform

# c-interface
rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

if platform.uname()[0] == "Windows":
    name = "libsolve_ivp.dll"
elif platform.uname()[0] == "Linux":
    name = "libsolve_ivp.so"
else:
    name = "libsolve_ivp.dylib"
libdop853 = ct.CDLL(rootdir+name)
dop853_solve_ivp_wrapper = libdop853.dop853_solve_ivp_wrapper

dop853_solve_ivp_wrapper.argtypes = [
    ct.c_void_p, # rhs_fcn_c
    ct.c_void_p, # t_span
    ct.c_int,    # neq
    ct.c_void_p, # y0
    ct.c_int,    # nt
    ct.c_void_p, # t_eval
    ct.c_int,    # n_events
    ct.c_void_p, # events_fcn
    ct.c_void_p, # data_ptr
    ct.c_double, # first_step
    ct.c_double, # max_step
    ct.c_double, # rtol
    ct.c_double, # atol
    ct.c_void_p, # len_message
    ct.c_void_p, # nt_out
    ct.c_void_p, # n_events_found
    ct.c_void_p, # d_p
]
dop853_solve_ivp_wrapper.restype = None

dop853_solve_ivp_result = libdop853.dop853_solve_ivp_result

dop853_solve_ivp_result.argtypes = [
    ct.c_void_p,
    ct.c_int,
    ct.c_int,
    ct.c_int,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
]
dop853_solve_ivp_result.restype = None


# ODE result
spec = [
    ('message', types.unicode_type),
    ('nfev', types.int64),
    ('success', types.boolean),
    ('t', types.double[:]),
    ('y', types.double[:,:]),
    ('ind_event', types.int64),
    ('t_event', types.double),
    ('y_event', types.double[:]),
]

@nb.experimental.jitclass(spec)
class ODEResult():
    def __init__(self, message, nfev, success, t, y, ind_event, t_event, y_event):
        self.message = message
        self.nfev = nfev
        self.success = success
        self.t = t
        self.y = y
        self.ind_event = ind_event
        self.t_event = t_event
        self.y_event = y_event

# User-facing interface
@nb.njit()
def solve_ivp(funcptr, t_span, y0, t_eval = np.array([],np.double), method = 'DOP853', \
          n_events = 0, event_fcn = 0, data = 0, \
          first_step = -1.0, max_step = -1.0, rtol = 1.0e-3, atol = 1.0e-6, min_step = 0.0):
    return _solve_ivp(funcptr, t_span, y0, t_eval, method, \
                    n_events, event_fcn, data, \
                    first_step, max_step, rtol, atol, min_step)

# Function with typing signature
signature = ODEResult.class_type.instance_type(
    types.int64, 
    types.Array(types.double, 1, 'C', readonly=True),
    types.Array(types.double, 1, 'C', readonly=True), 
    types.Array(types.double, 1, 'C', readonly=True),
    types.unicode_type,
    types.int32,
    types.int64,
    types.int64,
    types.double,
    types.double,
    types.double,
    types.double,
    types.double
)

@nb.njit(signature)
def _solve_ivp(funcptr, t_span, y0, t_eval, method, \
          n_events, event_fcn, data, \
          first_step, max_step, rtol, atol, min_step):

    if t_span.size != 2:
        raise ValueError('"t_span" must be length 2')

    neq = y0.size
    nt = t_eval.size

    len_message = np.empty((),np.int32)
    nt_out = np.empty((),np.int32)
    n_events_found = np.empty((),np.int32)
    d_p = np.empty((),np.int64)

    dop853_solve_ivp_wrapper(
        funcptr,
        t_span.ctypes.data,
        neq,
        y0.ctypes.data,
        nt,
        t_eval.ctypes.data,
        n_events,
        event_fcn,
        data,
        first_step,
        max_step,
        rtol,
        atol,
        len_message.ctypes.data,
        nt_out.ctypes.data,
        n_events_found.ctypes.data,
        d_p.ctypes.data
    )

    message_ = np.empty((len_message.item()+1),dtype=np.uint8)
    success = np.empty((),np.bool_)
    nfev = np.empty((),np.int32)
    t = np.empty((nt_out.item(),),np.double)
    y = np.empty((nt_out.item(),neq,),np.double)

    ind_event = np.ones((),np.int32)
    t_event = np.zeros((),np.double)
    y_event = np.zeros((neq,),np.double)

    dop853_solve_ivp_result(
        d_p.item(),
        neq,
        len_message.item(),
        nt_out.item(),
        message_.ctypes.data,
        success.ctypes.data,
        nfev.ctypes.data,
        t.ctypes.data,
        y.ctypes.data,
        ind_event.ctypes.data,
        t_event.ctypes.data,
        y_event.ctypes.data,
    )

    message = ''.join([chr(z) for z in message_])
    return ODEResult(message, nfev.item(), success, t, y, ind_event, t_event, y_event)



