c
c
c     Purpose: Fortran driver of genbtf, generate block triangular
c     form.  A structurally nonsingular matrix is permuted 
c
c     (A(rowperm(i), (colperm(j)), 1 <= i,j, <=n) 
c
c     to block upper triangular with sqcmpn blocks.
c     If matsng != 0, then the matrix is structurally singular
c     (for example has a zero row), and has a complex 
c     block triangular form.  
c


c
c      rcmstr ccmstr
c
c     If matsng != 0, then the matrix is structurally singular
c     (for example has a zero row), and has a complex 
c     block triangular form.  
c
c
c     test input
      parameter (n=3, nnz= 3)
c     and C-style indexing.
c
c
c
c
c
c     Start with an n by n matrix in sparse row format
c     ja, ia: column indices and row pointers
c
c     integer input scalars
c      integer msglvl = 0, output = 6
      integer msglvl, output
c
c     integer input arrays
      integer ia(n+1) , ja(nnz), iat(n+1), jat(nnz)
c
c     local work space
      integer w(10*n), i, j
c
c     integer output scalars
c 
      integer matsng
c     horizontal block:  rows, columns, connected components 
      integer nhrows, nhcols, hrzcmp
c     square block:  rows=columns, connected components 
      integer nsrows, sqcmpn
c     vertical block:  rows, columns, connected components 
      integer nvrows, nvcols, vrtcmp
c
c     integer output arrays
c     rowperm: row permutation, 
c     cotn: column permutation, old to new
      integer colperm(n), rowperm(n), rcmstr(n+1), ccmstr(n+1)
      matsng = 0
      msglvl = 0
      output = 6
c
c     More test input
      ia(1) = 0
      ia(2) = 1
      ia(3) = 2
      ia(4) = 3 
c
      ja(1) = 1
      ja(2) = 2 
      ja(3) = 0
c
c
c     Convert from C indexing to Fortran
c      if( nnz != ia(n+1) )then
c         stop         I can not remember Fortran syntax.
c      endif
      do 100 i=1,n+1
        ia(i) = ia(i) + 1
 100  continue
      do 101 i=1,nnz
        ja(i) = ja(i) + 1
 101  continue
      call mattrans(n,n,ja,ia,jat,iat)
c
c
      print*,'Input (row, column)'
      do 200 i=1,n
        do 201 j= ia(i),ia(i+1)-1
          print*,'    ',i,ja(j)
 201  continue
 200  continue
c
c
      call genbtf( n, n, 
     $             iat , jat,     ia,     ja, w     , 
     $             rowperm  , colperm  , nhrows,
     $             nhcols, hrzcmp, nsrows, sqcmpn, nvrows,
     $             nvcols, vrtcmp, rcmstr, ccmstr, msglvl, output )
c
c
      if( nhrows .gt. 0) then
        print*, "horizontal block:", nhrows, nhcols, hrzcmp
      endif
      print*, sqcmpn, "  blocks"
      if( nvrows .gt. 0) then
        print*, "vertical block:", nvrows, nvcols, vrtcmp
      endif
      matsng = nhrows + nvrows + hrzcmp + vrtcmp
      if( matsng .eq. 0) then
        do 401 i=1, sqcmpn
          print*,'    ', rcmstr(hrzcmp+i), ccmstr(hrzcmp+i)
 401    continue
      else
        print*, 'Structurally singular matrix'
      endif
c
      print*,'Permuted (row, column)'
      do 300 i=1,n
        k = rowperm(i)
        do 301 j= ia(k),ia(k+1)-1
          print*,'    ', i,ja(colperm(j))
 301  continue
 300  continue
c
c  rowperm --, colperm --
c  rcmstr --, ccmstr --
c
c
c
c
      end
