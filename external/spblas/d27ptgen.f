*--------------------------------------------------------------------
*\Documentation
*
*\Name: S27PGEN
*
*\Description:
*    S27PGEN generates a matrix stored in (block) compressed row format where
*    the matrix pattern is that of a 27-point finite difference
*    operator on a rectangular parallelepiped.  The values are random
*    between zero and one for the off-diagonal elements.  The main
*    diagonal is defined so that the matrix is diagonally-dominant.
*    The matrix can be generated in symmetric or non-symmetric form.
*    If a symmetric matrix is formed, only the lower triangular part
*    of the matrix is generated.
*
*
*\Arguments
*
*   nnzmx     integer (input)
*             On entry, nnzmx is the number of elements available in
*             ia,ja and nnzmx*nb*nb is the size of a.
*             nnzmx is unchanged on exit.
*
*    nx       integer (input)
*             On entry, nx is the number of gridpoints in x direction 
*             when generating PDE-type matrix.
*
*    ny       integer (input)
*             On entry, ny is the number of gridpoints in y direction 
*             when generating PDE-type matrix.
*
*    nz       integer (input)
*             On entry, nz is the number of gridpoints in z direction 
*             when generating PDE-type matrix.
*
*   nb        integer (input)
*             On entry, nb is the block entry size.  Each matrix entry
*             is of size (nb, nb).
*
*   factor    real (input)
*             On entry, factor is a real number used to define the 
*             strictly upper triangle of A as factor times the transpose 
*             of the strictly lower triangle of A.  This value is used
*             only if isym = 0.
*             factor is unchanged on exit.
*       
*   isym      integer (input)
*             On entry, isym indicates if the matrices generated should 
*             be symmetric or non-symmetric.
*              0 - nonsymmetric
*              1 - symmetric
*             isym is unchanged on exit.
*
*    title    character*80 (output)
*             On exit, title is an 80-character string describing 
*             the current matrix.  The first character is either 
*             an 'S' or a 'U' indicating whether or not the matrix 
*             is symmetric or unsymmetric, resp.
*
*    n        integer (output)
*             On exit, n is the requested block dim of the matrix to be
*             generated.  It is defined to be n = nx * ny * nz.
*
*    ia       integer array (output)
*             On entry, ia is an nnzmx-length integer array.  On exit,
*             the first nnz element are defined such that ia(i) is
*             the row position of the element a(i) in the matrix.
*
*    a        real array (output)
*             On entry, a is an (nb,nb,nnzmx) dim scalar array.  On exit,
*             the  (nb,nb,nnz) elements are the nonzero (block) elements of the
*             generated matrix.
*
*    ia       integer array (output)
*             On entry, ia is an nnzmx-length integer array.  On exit,
*             the first nnz element are defined such that ia(i) is
*             the row position of the element a(i) in the matrix.
*
*    ja       integer array (output)
*             On entry, ja is an nnzmx-length integer array.  On exit,
*             the first nnz element are defined such that ja(i) is
*             the column position of the element a(i) in the matrix.
*
*   nnz       integer (output)
*             On exit, nnz is the number of (block) nonzeros in the generated
*             matrix or, if isym = 1, it is the number of elements in
*             lower triangle of the matrix.
*
*
*   ierr      integer (output)
*             On exit, ierr contains the error flag.
*              0 - No error.
*            - 1 - Ran out of workspace.  Rerun with nnzmx larger.
*\Remarks:
*     None.
*\Examples:
*
*                      Data Structure
*                      --------------
*  The format used to generate the matrix is a sparse coordinate
*format that allows for representation of a large, sparse,
*unpatterned matrix.  The examples below illustrate how a matrix is
*represented in this format.
*
*          Symmetric case
*          --------------
*
*          Suppose the matrix A is symmetric with 
*
*                  | 11   0   0  14   0 |
*                  |                    |
*                  |  0  22  23   0  25 |
*                  |                    |
*              A = |  0  23  33  34   0 |.
*                  |                    |
*                  | 14   0  34  44   0 |
*                  |                    |
*                  |  0  25   0   0  55 |
*
*A is represented by the three arrays a, ia, ja.  In the
*symmetric case only the lower triangle of A is stored.  a contains
*the non-zero values of A not necessarily ordered in any way.
*ia contains the row indices and ja contains the column indices.
*
*In the above case, 
*
*         a = | 11  14  22  23  25  33  34  44  55 |,
*
*        ia = |  1   4   2   3   5   3   4   4   5 |
*
*and
*        ja = |  1   1   2   2   2   3   3   4   5 |
*
*Note that the order of the elements in a is not important as long
*as it corresponds to the order in ia and ja.
*
*          Non-symmetric case
*          ------------------
*
*          Suppose the matrix A is non-symmetric with 
*
*                  | 11   0   0  14   0 |
*                  |                    |
*                  |  0  22  23   0  25 |
*                  |                    |
*              A = |  0  32  33  34   0 |.
*                  |                    |
*                  | 41   0  43  44   0 |
*                  |                    |
*                  |  0  52   0   0  55 |
*
*In this case all the non-zero element of A are stored in the same
*manner.
*Thus,
*
*        a = | 11 41 22 32 52 23 33 43 14 34 44 25 55 |,
*
*       ia = |  1  4  2  3  5  2  3  4  1  3  4  2  5 |
*
*and
*
*       ja = |  1  1  2  2  2  3  3  3  4  4  4  5  5 |.
*
*\Enddoc
*
*--------------------------------------------------------------------
*
*\Lib
*
*\Local variables:
*
*    ichar0   integer
*             ASCII value of the character '0'.
*    charn    character*5
*             5-digit character representation of the matrix size 
*             (the argument n).
*    k, k1    integer
*             Implicit indices used to keep track of the current non-
*             zero element.
*    krow     Integer
*             Current matrix row.
*    i,j,l    Integer
*             Loop indices over the x,y and z direction, resp.
*
*
*\Author:
*    Mike Heroux
*
*\References:
*    None.
*
*\Keywords:
*    sparse matrices; PDE matrices
*
*\Routines called:
*    Intrinsics - mod, char, ranf()
*
*     \Implementation details:
*     None.
*     
*     \Performance data:
*     None.
*     
*     \Endlib
*     
*--------------------------------------------------------------------
      subroutine d27ptgen 
     &     ( nnzmx, nx, ny, nz, nb, factor, isym, 
     &     title, n, a, ia, ja,
     &     nnz, ierr )
*     
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer idump
      parameter(idump = 0)
      integer nnzmx, nnz, nx, ny, nz, nb,
     &     ia(*), ja(nnzmx), isym, ierr
      double precision a(nb,nb,nnzmx), factor
      character*80 title
      character*5 charn
*     
*     ----------------------------
*     Specifications for local vars
*     ----------------------------
      integer k, k1, krow, n, ichar0, itmp, i
      integer itmp1, l, j, iblk, jblk, ii
*     
*     ----------------------------
*     Specifications for commons
*     ----------------------------
      common k, k1, krow
*     
*     --------------------------
*     First executable statement
*     --------------------------
      ierr = 0
c.....define system dimension
      n = nx * ny * nz
c     
c.....convert n to character string if n < 10 000 000
c     
c.....get ASCII code for 0
      ichar0 = ichar('0')
      if (n.lt.100 000) then
         itmp = 1
         do i=5,1,-1
            itmp1 = mod((n/itmp),10)
            charn(i:i) = char(ichar0+itmp1)
            itmp = itmp*10
         end do
      elseif (n.lt.10 000 000) then
         itmp = 1 000
         do i=4,1,-1
            itmp1 = mod((n/itmp),10)
            charn(i:i) = char(ichar0+itmp1)
            itmp = itmp*10
         end do
         charn(5:5) = 'K'
      else
         charn = '*****'
      endif
c     
c.....check if we do symmetric or non-symmetric matrix
      if (isym.eq.0) then
c     
c.....***** Non-symmetric case *****
c     
c.....define title
         title = 'UNSYMMETRIC 27 POINT MATRIX N = ' // charn
         title(73:80) = 'U27' //charn
c     
c.....initialize counters
         k = 1
         k1 = 0
         krow = 0
c     
c.....loop through domain to get natural ordering
         do 5 l=1,nz
            do 10 j=1,ny
               do 20 i=1,nx
                  k = k + k1
                  k1 = 0
                  krow = krow + 1
                  ia(krow) = k
                  if (k+27.gt.nnzmx) then
c...........ran out of space so set error flag and return
                     ierr = - 1
                     goto 9999
                  endif
c     
c.........reach to left
                  if (i.gt.1) then
                     call update(a,ia,ja,nb,krow - 1)
                     if (j.gt.1) then
                        call update(a,ia,ja,nb,krow - nx - 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny-nx-1)
                        endif
                        if (l.lt.nz) then
                           call update(a,ia,ja,nb,krow+nx*ny-nx-1)
                        endif
                     endif
                     if (j.lt.ny) then
                        call update(a,ia,ja,nb,krow + nx - 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny+nx-1)
                        endif
                        if (l.lt.nz) then
                           call update(a,ia,ja,nb,krow+nx*ny+nx-1)
                        endif
                     endif
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny - 1)
                     endif
                     if (l.lt.nz) then
                        call update(a,ia,ja,nb,krow + nx * ny - 1)
                     endif
                  endif
c     
c.........reach below
                  if (j.gt.1) then
                     call update(a,ia,ja,nb,krow - nx)
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny - nx)
                     endif
                     if (l.lt.nz) then
                        call update(a,ia,ja,nb,krow + nx * ny - nx)
                     endif
                  endif
c     
c.........reach back
                  if (l.gt.1) then
                     call update(a,ia,ja,nb,krow - nx * ny)
                  endif
c     
c.........reach to right
                  if (i.lt.nx) then
                     call update(a,ia,ja,nb,krow + 1)
                     if (j.gt.1) then
                        call update(a,ia,ja,nb,krow - nx + 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny-nx+1)
                        endif
                        if (l.lt.nz) then
                           call update(a,ia,ja,nb,krow+nx*ny-nx+1)
                        endif
                     endif 
                     if (j.lt.ny) then
                        call update(a,ia,ja,nb,krow + nx + 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny+nx+1)
                        endif
                        if (l.lt.nz) then
                           call update(a,ia,ja,nb,krow+nx*ny+nx+1)
                        endif
                     endif 
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny + 1)
                     endif
                     if (l.lt.nz) then
                        call update(a,ia,ja,nb,krow + nx * ny + 1)
                     endif
                  endif
c     
c.........reach above
                  if (j.lt.ny) then
                     call update(a,ia,ja,nb,krow + nx)
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny + nx)
                     endif
                     if (l.lt.nz) then
                        call update(a,ia,ja,nb,krow + nx * ny + nx)
                     endif
                  endif
c     
c.........reach forward
                  if (l.lt.nz) then
                     call update(a,ia,ja,nb,krow + nx * ny)
                  endif
c     
c.........fill in diagonal term
                  do jblk=1,nb
                     do iblk=1,nb
                        a(iblk,jblk,k+k1) = 0.0
                        do 30 ii = 0,k1-1
                           a(iblk,jblk,k+k1) = a(iblk,jblk,k+k1) - 
     &                          a(iblk,jblk,k+ii)
 30                     continue
                     end do
                  end do
c     
c.........make diagonally dominant
                  do iblk=1,nb
                     a(iblk,iblk,k+k1) = 2.0*a(iblk,iblk,k+k1)
                  end do
                  ja(k+k1) = krow
                  k1 = k1 + 1
 20            continue
 10         continue
 5       continue
c     
c.....set number of non-zero elements
         nnz = k + k1 - 1
         ia(krow+1) = nnz+1
c     
      else
c     
c.....symmetric case
c     
c.....define title
         title = 'SYMMETRIC 27 POINT MATRIX N = ' // charn
         title(73:80) = 'S27' //charn
c     
c     
c.....initialize counters
         k = 1
         k1 = 0
         krow = 0
c     
c.....loop through domain to get natural ordering
         do 105 l=1,nz
            do 110 j=1,ny
               do 120 i=1,nx
                  k = k + k1
                  k1 = 0
                  krow = krow + 1
                  ia(krow) = k
                  if (k+14.gt.nnzmx) then
c...........ran out of space so set error flag and return
                     ierr = - 1
                     goto 9999
                  endif
c     
c.........reach to left
                  if (i.gt.1) then
                     call update(a,ia,ja,nb,krow - 1)
                     if (j.gt.1) then
                        call update(a,ia,ja,nb,krow - nx - 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny-nx-1)
                        endif
                     endif
                     if (j.lt.ny) then
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny+nx-1)
                        endif
                     endif
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny - 1)
                     endif
                  endif
c     
c.........reach below
                  if (j.gt.1) then
                     call update(a,ia,ja,nb,krow - nx)
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny - nx)
                     endif
                  endif
c     
c.........reach back
                  if (l.gt.1) then
                     call update(a,ia,ja,nb,krow - nx * ny)
                  endif
c     
c.........reach to right
                  if (i.lt.nx) then
                     if (j.gt.1) then
                        call update(a,ia,ja,nb,krow - nx + 1)
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow-nx*ny-nx+1)
                        endif
                     endif 
                     if (j.lt.ny) then
                        if (l.gt.1) then
                           call update(a,ia,ja,nb,krow - nx*ny+nx+1)
                        endif
                     endif 
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny + 1)
                     endif
                  endif
c     
c.........reach above
                  if (j.lt.ny) then
                     if (l.gt.1) then
                        call update(a,ia,ja,nb,krow - nx * ny + nx)
                     endif
                  endif
c     
c.........fill in diagonal term
                  do jblk=1,nb
                     do iblk=1,nb
                        a(iblk,jblk,k+k1) = 0.0
                        do 130 ii = 0,k1-1
                           a(iblk,jblk,k+k1) = a(iblk,jblk,k+k1) - 
     &                          a(iblk,jblk,k+ii)
 130                    continue
                     end do
                  end do
c     c     
c.........make diagonally dominant
                  do iblk=1,nb
                     a(iblk,iblk,k+k1) = 14.0 + 1.001*a(iblk,iblk,k+k1)
                  end do
                  ja(k+k1) = krow
                  k1 = k1 + 1
 120           continue
 110        continue
 105     continue
c     
c.....set number of non-zero elements
         nnz = k + k1 - 1
         ia(krow+1) = nnz+1
c     
      endif
c     
 9999 continue
*
*     This code is useful for debugging.  It dumps out the matrix in
*     coordinate format: i, j, A(i,j) which can be read by Matlab.
*     This is turned on by setting the parameter idump = 1 at top of routine.
*
      if (idump.eq.1) then
         open(11,file='temp.dat')
         do i = 1, n
            do j = ia(i),ia(i+1)-1
               do iblk=1,nb
                  do jblk=1,nb
                     write(11,98)
     &               (i-1)*nb+iblk,(ja(j)-1)*nb+jblk,a(iblk,jblk,j)
                  end do
               end do
            end do
         end do
         close(11)
      endif
 98   format(i10,1x, i10, 1x, e26.16)
      return
      end
*     
*     Helper function for filling off-diagonal terms
*     
      subroutine update(a,ia,ja,nb,icol)
      implicit none
      integer ia(*), ja(*), k, k1, krow, nb
      integer icol, jblk, iblk
      double precision a(nb,nb,*)
      double precision ranf_local
      common k, k1, krow
*     
      do jblk=1,nb
         do iblk=1,nb
            a(iblk,jblk,k+k1) = -ranf_local()
         end do 
      end do
      ja(k+k1) = icol
      k1 = k1 + 1
*      print*,'k,k1,icol = ',k,k1,icol
*      if (krow+icol.lt.0) then
*      print*,'krow+icol=',krow+icol
*      endif
      return
      end
