import ctypes
from pathlib import Path
from pkg_resources import resource_filename, resource_listdir

import numpy as np


compiled_libraries = {
    filename.split(".")[0]: Path(resource_filename("fundem", filename))
    for filename in resource_listdir("fundem", ".")
    if filename.startswith("_") and filename.endswith(".so")
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
    ctypes.int,
    ctypes.c_size_t
]


def first_moment_survival(mx, ax, nx):
    assert mx.shape == ax.shape
    assert mx.shape[-1] == nx.shape[-1], \
        f"Populations {mx.shape} must end in age groups {nx.shape}."
    population_cnt = 1 if len(mx.shape) is 1 else np.prod(mx.shape[:-1])

    survival = np.require(np.zeros_like(mx), dtype=np.float64,
                          requirements=["A", "C", "W", "O"])
    _lifetable.first_moment_survival(
        np.ascontiguousarray(mx).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(ax).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(nx).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        survival.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        nx.size,
        population_cnt
    )
    return survival


_lifetable.first_moment_population.restype = None
_lifetable.first_moment_population.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_size_t
]


def first_moment_population(mx, ax, nx):
    assert mx.shape == ax.shape
    assert mx.shape[-1] == nx.shape[-1], \
        f"Populations {mx.shape} must end in age groups {nx.shape}."
    population_cnt = 1 if len(mx.shape) is 1 else np.prod(mx.shape[:-1])

    lx = np.require(np.zeros_like(mx), dtype=np.float64,
                    requirements=["A", "C", "W", "O"])
    dx = np.require(np.zeros_like(mx), dtype=np.float64,
                    requirements=["A", "C", "W", "O"])
    _lifetable.first_moment_population(
        np.ascontiguousarray(mx).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(ax).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(nx).ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)),
        lx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        dx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        nx.size,
        population_cnt
    )
    return lx, dx
