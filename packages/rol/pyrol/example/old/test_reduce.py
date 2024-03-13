from PyROL.PyROL import ROL

op = ROL.Elementwise.ReductionMin_double_t()

v = 10.
w = 20.
minVal = min(v, w)
op.reduce(v, w)
assert abs(w-minVal) < 1e-7, (v, w, minVal)
