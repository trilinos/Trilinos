from pylab import *

c = 'bgr'
y = array([0, 1, 0])

# PE 0
i = 0
x0 = array([0, 1])
x1 = array([0, 1, 2])
x2 = array([1, 2])
x3 = array([2, 3])
plot(x0, y[1:], c[i], x1, y, c[i], x2, y[:2], c[i], x3, y[:2], c[i])

i = 1
x0 = array([2, 3])
x1 = array([3, 4])
x2 = array([3, 4, 5])
x3 = array([4, 5])
x4 = array([5, 6])
plot(x0, y[1:], c[i], x1, y[1:], c[i], x2, y, c[i], x3, y[:2], c[i], x4, y[:2], c[i])

i = 2
x0 = array([5, 6])
x1 = array([6, 7])
x2 = array([6, 7, 8])
x3 = array([7, 8])
plot(x0, y[1:], c[i], x1, y[1:], c[i], x2, y, c[i], x3, y[:2], c[i])

yticks([], [])
savefig('../Images/export_hatfun.pdf', transparent=True)
