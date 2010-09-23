M1 = [ 1 0 2 ; 2 5 0 ; 0 3 1 ];
M2 = [ 5 0 6 ; 2 6 0 ; 0 3 7 ];

A = [ 1.0*M2 2.0*M1 ; 3.0*M1 4.0*M2 ];
Qp = [ 9 8 7 ; 1 2 3 ; 4 5 6];
iQp = diag(1./diag(Qp));

Lap = M2*M1;
invA00 = inv(M2)
invS = inv(Lap)*M1*iQp
