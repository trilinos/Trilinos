A = [  1  2  1 -1
       2  1 -3  1
       1 -3  1  2
      -1  1  2  1 ];

% test_result
v = [ 0  1  1  3]';
(inv(A)*v)'

v = [ -2 4 7 9]'; 
(inv(A)*v)'

v = [ 1 0 0 -5]';
(inv(A)*v)'

v = [ 4 -4 6 12]';
(inv(A)*v)'

% test_alphabeta
alpha = 0.5;
beta = -1.9;

u = [6.2 -2.4 9.7 0.04]';

v = [ 0  1  1  3]';
(alpha * inv(A)*v + beta*u)'

v = [ -2 4 7 9]'; 
(alpha * inv(A)*v + beta*u)'

v = [ 1 0 0 -5]';
(alpha * inv(A)*v + beta*u)'

v = [ 4 -4 6 12]';
(alpha * inv(A)*v + beta*u)'
