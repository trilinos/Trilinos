      program f77_main
      integer maxn, maxnz
      parameter ( maxn = 10000, maxnz = 100000 )
      integer rowind(maxnz), colptr(maxn)
      real*8  values(maxnz), b(maxn)
      integer n, nnz, nrhs, ldb, info
*
      call hbcode1(n, n, nnz, values, rowind, colptr)
*
      nrhs = 1
      ldb = n
      do i = 1, n
         b(i) = 1
      enddo
*
      call c_bridge_dgssv( n, nnz, nrhs, values, rowind, colptr, 
     $                     b, ldb, info )
*
      if (info .eq. 0) then
          write (*,*) (b(i), i=1, 10)
      else
         write(*,*) 'INFO from c_bridge_dgssv = ', info
      endif

      stop
      end


