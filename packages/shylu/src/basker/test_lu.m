
index = UFget ;

% get all square matrices of size less than 1000
f = find (index.nrows == index.ncols & index.nrows < 1001 & index.isReal & ~index.isBinary);

[ignore,i] = sort (index.nrows (f)) ;
f = f (i) ;

skiplist = [ 905 ] ;

fid = fopen('test.data', 'w');

for i = f
  fprintf (fid, '%d: %s/%s: ', i, index.Group {i}, index.Name{i}) ;
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

    % Set BTF to off for fair comparison.
    opts.btf = 0;

    tic
    [L, U, P, Q] = lu (A) ;
    matlab_t = toc;

    B = P*A*Q ;

    if (nnz (diag (U)) == n)
        % LU factorization worked
        % do your solve.  It should not return x sorted,
        % but only in a correct topological order

        tic
        [LU, info] = klu(B, opts) ;
        klu_t = toc;

        tic
        [L1, U1, P1] = basker(B) ;
        basker_t = toc;

        nn = norm (L1*U1 - P1*B, 1)/norm(B, 1);
fprintf (fid,'                        your   resid: %g klu_time:%g basker_time:%g speedup: %g\n', nn, klu_t, basker_t, klu_t/basker_t) ;
        % compare with MATLAB
        %xgood = L\b ;
        %norm (full (x-xgood)) ;
        %fprintf ('          your   resid: %g\n', normest (L1*U1-B)) ;
        %input(' Hit enter to continue') ;
        %fprintf ('          MATLAB resid: %g\n', norm (L*xgood-b)) ;
    else
      fprintf (fid, '	    skip (zero pivot encountered)\n') ;
    end

    clear L  A U P Q B L1 U1 LU
    %pack

end

    fclose(fid);
