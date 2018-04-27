import numpy as np
cimport numpy as np
import cython

cdef extern from "math.h":
    double abs(double t)

@cython.wraparound(False)
@cython.boundscheck(False)
def pairwise_distance(np.ndarray[np.double_t, ndim=1] r):
    cdef int i, j, c, size
    cdef np.ndarray[np.double_t, ndim=1] ans
    size = sum(range(1, r.shape[0]+1))
    ans = np.empty(size, dtype=r.dtype)
    c = -1
    for i in range(r.shape[0]):
        for j in range(i, r.shape[0]):
            c += 1
            ans[c] = abs(r[i] - r[j])
    return ans
