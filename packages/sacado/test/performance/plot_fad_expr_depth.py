from matplotlib.pylab import *
from numpy import *
from scipy.linalg import *

dfad = loadtxt('fad_expr_depth_dfad.txt')
elr = loadtxt('fad_expr_depth_elr_dfad.txt')
cache = loadtxt('fad_expr_depth_cache_dfad.txt')
elr_cache = loadtxt('fad_expr_depth_elr_cache_dfad.txt')

mult = [1, 2, 3, 4, 5, 10, 15, 20]
P = mult
i1 = 3
i2 = 9
n = 10

figure(1,figsize=(13,9));
clf();
subplot(2,1,1)
plot(P,dfad[i1,1:], 'k-s',
     P,dfad[i2,1:], 'k--s')
legend(('Standard N = 5', 'Standard N = 50'), loc='upper left')
ylabel('Scaled Run Time')
title('Mult')
subplot(2,1,2)
plot(P,elr[i1,1:], 'k-s',
     P,cache[i1,1:], 'k--*',
     P,elr_cache[i1,1:], 'k:D',
     P,elr[i2,1:], 'k-p',
     P,cache[i2,1:], 'k--h', 
     P,elr_cache[i2,1:], 'k:>')
legend(('ELR N = 5', 'Cache N = 5', 'Cache ELR N = 5', 'ELR N = 50', 'Cache N = 50', 'Cache ELR N = 50'), loc='upper left')
xlabel('Expression Size P')
ylabel('Scaled Run Time')
savefig('mult.pdf')
savefig('mult.eps')

figure(2,figsize=(13,9));
clf();
subplot(2,1,1)
plot(P,dfad[i1+n,1:], 'k-s',
     P,dfad[i2+n,1:], 'k--s')
legend(('Standard N = 5', 'Standard N = 50'), loc='upper left')
ylabel('Scaled Run Time')
title('Add')
subplot(2,1,2)
plot(P,elr[i1+n,1:], 'k-s',
     P,cache[i1+n,1:], 'k--*',
     P,elr_cache[i1+n,1:], 'k:D',
     P,elr[i2+n,1:], 'k-p',
     P,cache[i2+n,1:], 'k--h', 
     P,elr_cache[i2+n,1:], 'k:>')
legend(('ELR N = 5', 'Cache N = 5', 'Cache ELR N = 5', 'ELR N = 50', 'Cache N = 50', 'Cache ELR N = 50'), loc='upper left')
xlabel('Expression Size P')
ylabel('Scaled Run Time')
savefig('add.pdf')
savefig('add.eps')

figure(3,figsize=(13,9));
clf();
subplot(2,1,1)
plot(P,dfad[i1+2*n,1:], 'k-s',
     P,dfad[i2+2*n,1:], 'k--s')
legend(('Standard N = 5', 'Standard N = 50'), loc='upper left')
ylabel('Scaled Run Time')
title('Nested')
subplot(2,1,2)
plot(P,elr[i1+2*n,1:], 'k-s',
     P,cache[i1+2*n,1:], 'k--*',
     P,elr_cache[i1+2*n,1:], 'k:D',
     P,elr[i2+2*n,1:], 'k-p',
     P,cache[i2+2*n,1:], 'k--h', 
     P,elr_cache[i2+2*n,1:], 'k:>')
legend(('ELR N = 5', 'Cache N = 5', 'Cache ELR N = 5', 'ELR N = 50', 'Cache N = 50', 'Cache ELR N = 50'), loc='upper right',bbox_to_anchor=(1.0,1.8))
xlabel('Expression Size P')
ylabel('Scaled Run Time')
savefig('nested.pdf')
savefig('nested.eps')

