import ctypes
from pathlib import Path
from pkg_resources import resource_filename, resource_listdir

import numpy as np


compiled_libraries = {
    filename.split(".")[0]: Path(resource_filename("fundem", filename))
    for filename in resource_listdir("fundem", ".")
    if filename.startswith("_")
}

_lifetable = np.ctypeslib.load_library(
    compiled_libraries["_lifetable"].name,
    compiled_libraries["_lifetable"].parent)

_lifetable.first_moment_survival.restype = None
_lifetable.first_moment_survival.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
]
def first_moment_survival(mx, ax, nx):
    survival = np.require(np.zeros_like(mx), dtype=np.float64,
                          requirements=["A", "C", "W", "O"])
    _lifetable.first_moment_survival(
        np.ascontiguousarray(mx).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(ax).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(nx).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        survival.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        mx.size
    )
    return survival
