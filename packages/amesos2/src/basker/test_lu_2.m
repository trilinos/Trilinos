
index = UFget ;

% get all square matrices of size less than 1000
f = find (index.nrows == index.ncols & index.nrows < 1001 & index.isReal & ~index.isBinary);

[ignore,i] = sort (index.nrows (f)) ;
f = f (i) ;

skiplist = [ 905 ] ;

fid = fopen('test.data', 'w');
fid = 1;
for i = f
  fprintf (fid, '%d: %s/%s: ', i, index.Group {i}, index.Name{i}) ;
	
   %pause;
  %input(' Hit enter to continue') ;
    Problem = UFget (i) ;
    A = Problem.A ;
    clear Problem
    n = size (A,1) ;
    [m , n] = size(A) ;

    if (any (skiplist == i))
        fprintf ('	    skip\n') ;
        continue
    end

    
    %get rhs%
    x = ones(n,1);
    y = A*x;
     
    tic;
    [L1, U1, P1, out] = basker(A, nnz(A), nnz(A), y) ;
    basker_t = toc;
    nn = norm (x - out, 1)/norm(x, 1);
    fprintf (fid,'  your   resid: %g ', nn) ;
        % compare with MATLAB
        %xgood = L\b ;
        %norm (full (x-xgood)) ;
        %fprintf ('          your   resid: %g\n', normest (L1*U1-B)) ;
        %input(' Hit enter to continue') ;
        %fprintf ('          MATLAB resid: %g\n', norm (L*xgood-b)) ;
   
      fprintf (fid, '	    skip (zero pivot encountered)\n') ;
    

    clear L  A U P Q B L1 U1 LU
    %pack

end

    fclose(fid);
