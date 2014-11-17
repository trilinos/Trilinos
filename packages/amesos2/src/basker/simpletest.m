%Pr = UFget(905) ;
Pr = UFget(449) ;
%Pr = UFget(23) ;
A = Pr.A ;
%[ L, U, P, Q] = lu(A) ;
%B=P*A*Q ;
%fprintf(' Checking max memory \n ' ) ;
%[L, U] = basker(B) ;
%norm=normest(L*U-B) 
%fprintf(' Test 1 Done \n ' ) ;

%A = delsq(numgrid('S', 10));
[L,U,Q] = basker(A, nnz(A), nnz(A));
