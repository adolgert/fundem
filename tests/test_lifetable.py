import numpy as np

from fundem import lifetable


def test_first_moment_survival():
    N = 23
    mx = np.full((N,), 0.1, dtype=np.float)
    ax = np.full((N,), 2.5, dtype=np.float)
    nx = np.full((N,), 5, dtype=np.float)
    survival = lifetable.first_moment_survival(mx, ax, nx)
    assert len(survival) == N
    assert np.all(survival < 1.0)
    assert np.all(survival > 0)
