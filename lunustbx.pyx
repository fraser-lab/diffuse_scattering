import numpy as np
cimport numpy as np
cimport cython
from cython.parallel import prange, parallel
from libc.math cimport M_PI, sqrt, atan, cos, sin, acos


DTYPE = np.int64
# ctypedef np.int_t DTYPE_t

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
def thrsh_filt(long[:,::1] inp, long minval, long maxval):
    cdef long n = inp.shape[0]
    cdef long m = inp.shape[1]
    cdef long x, y
    cdef long mask
    cdef long[:,::1] outp = np.empty_like(inp, dtype=DTYPE)
    mask = -32768
    with nogil:
        for y in prange(n):
            for x in range(m):
                if inp[y,x]==mask:
                    outp[y,x]=mask
                else:
                    if inp[y,x] < minval:
                        outp[y,x] = mask
                    elif inp[y,x] > maxval:
                        outp[y,x] = mask
                    else:
                        outp[y,x] = inp[y,x]
    return outp


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
def beamstop_filt(long[:,::1] inp, long xmin, long xmax, long ymin, long ymax):
    cdef long n = inp.shape[0]
    cdef long m = inp.shape[1]
    cdef long x, y
    cdef long mask
    cdef long[:,::1] outp = np.empty_like(inp, dtype=DTYPE)
    mask = -32768
    with nogil:
        for y in prange(n):
            for x in range(m):
                if inp[y,x]==mask:
                    outp[y,x]=mask
                else:
                    if x < xmin:
                        outp[y,x] = inp[y,x]
                    elif x > xmax:
                        outp[y,x] = inp[y,x]
                    else:
                        if y < ymin:
                            outp[y,x] = inp[y,x]
                        elif y > ymax:
                            outp[y,x] = inp[y,x]
                        else:
                            outp[y,x] = mask
    return outp


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
def edge_filt(long[:,::1] inp, long xmin, long xmax, long ymin, long ymax):
    cdef long n = inp.shape[0]
    cdef long m = inp.shape[1]
    cdef long x, y
    cdef long mask
    cdef long[:,::1] outp = np.empty_like(inp, dtype=DTYPE)
    mask = -32768

    with nogil:
        for y in prange(n):
            for x in range(m):
                if inp[y,x]==mask:
                    outp[y,x]=mask
                else:
                    if x < xmin:
                        outp[y,x] = mask
                    elif x > xmax:
                        outp[y,x] = mask
                    else:
                        if y < ymin:
                            outp[y,x] = mask
                        elif y > ymax:
                            outp[y,x] = mask
                        else:
                            outp[y,x] = inp[y,x]

    return outp


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
@cython.cdivision(True)
def polar_filt(long[:,::1] inp, long beam_x, long beam_y, double det_dist, double polar_fac, double polar_offset, double pixel_size):
    cdef long N = inp.shape[0]
    cdef long M = inp.shape[1]
    cdef long mask
    cdef long x, y
    cdef long r_y, r_x, radius_squared
    cdef double radius, distance_pixels, arctan_argument, cos_two_theta, sin_two_theta, cos_two_rho
    cdef long[:,::1] outp = np.empty_like(inp, dtype=DTYPE)
    cdef double POLARIZATION_CORRECTION_THRESHOLD = 0.01
    mask = -32768
    cdef double two_rho_offset, two_theta
    cdef double fact
    cdef long result
    
    two_rho_offset = 2.*M_PI/180.*polar_offset
    with nogil:
        for y in prange(N):
            r_y = (y - beam_y)
            for x in range(M):
                if inp[y,x] == mask:
                    outp[y,x] = mask
                else:
                    r_x = (x - beam_x)
                    radius_squared = ((r_x*r_x) + (r_y*r_y))
                    radius = sqrt(radius_squared)
                    distance_pixels = det_dist / pixel_size
                    arctan_argument = radius / distance_pixels
                    if (arctan_argument > POLARIZATION_CORRECTION_THRESHOLD):
                        two_theta = atan(arctan_argument)
                        cos_two_theta = cos(two_theta)
                        sin_two_theta = sin(two_theta)
                        cos_two_rho = cos(2*acos(r_x / radius) - two_rho_offset)
                        fact=(inp[y,x] * 2. / (1. + cos_two_theta*cos_two_theta -
                                    polar_fac*cos_two_rho * sin_two_theta*sin_two_theta))
                        result = (<long>fact)
                    outp[y,x]=result

    return outp


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
@cython.cdivision(True)
def solid_angle_filt(long[:,::1] inp, long beam_x, long beam_y, double det_dist, double tilt_x, double tilt_y, double pixel_size):
    cdef long N = inp.shape[0]
    cdef long M = inp.shape[1]
    cdef long mask
    cdef long x, y
    cdef double r_y, r_x
    cdef double radius_squared
    cdef long[:,::1] outp = np.empty_like(inp, dtype=DTYPE)
    cdef double correction_factor
    mask = -32768
    cdef double distance_squared
    cdef double fact
    cdef long result
    
    distance_squared = det_dist*det_dist
    
    with nogil:
        for y in prange(N):
            r_y = (y - beam_y)*pixel_size
            for x in range(M):
                if inp[y,x] == mask:
                    outp[y,x] = mask
                else:
                    r_x = (x - beam_x)*pixel_size
                    radius_squared = ((r_x*r_x) + (r_y*r_y))
                    correction_factor = (1. + radius_squared/distance_squared - 
                                         2*M_PI/180.*(r_x*tilt_y -r_y*tilt_x)/det_dist)
                    fact=(inp[y,x] * correction_factor * sqrt(correction_factor))
                    result = (<long>fact)
                    outp[y,x]=result

    return outp


# define a function pointer to a metric
# ctypedef long (*metric_ptr)(long[::1])
# ctypedef long (*fast_ptr)(long[:, ::1])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef inline long mode_finder(long[::1] X, long best=0) nogil:
    cdef long k
    cdef long mode
    cdef long i = X.shape[0]
    for k in range(i):
        if X[k] > best:
            best = X[k]
            mode = k
    return mode

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
def mode_filt(long[:, ::1] inp, long filt_n, long max_int):
    cdef:
        long N = inp.shape[0]
        long M = inp.shape[1]
        long x, y, xmin, xmax, ymin, ymax, xx, yy
        long mask = -32768
        long[::1] modc = np.empty(max_int, dtype=DTYPE)
        long[:, ::1] outp = np.empty_like(inp, dtype=DTYPE)
    
    for y in prange(N, nogil=True):
        for x in range(M):
            if inp[y,x] == mask:
                outp[y,x] = mask
            else:
                xmin = x - filt_n
                xmax = x + filt_n + 1
                ymin = y - filt_n
                ymax = y + filt_n + 1
                modc[:] = 0

                for yy in range(ymin,ymax):
                    if yy < 0:
                        pass
                    elif yy > (N-1):
                        pass
                    else:
                        for xx in range(xmin,xmax):
                            if xx < 0:
                                pass
                            elif xx > (M-1):
                                pass
                            else:
                                if inp[yy,xx] == mask:
                                    pass
                                else:
                                    modc[inp[yy,xx]]+= 1
                outp[y,x] = mode_finder(modc)
    
    
    return outp