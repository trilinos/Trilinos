
      subroutine initcos(n, x, mmax)
      implicit none 
      double precision x, PI
      integer n, mmax, k, j
      dimension x(0:mmax, 0:mmax)
      PI = 3.14159265358979d0
      do k=0, n
        do j=0, n
          x(k,j)=cos(PI*k*j/dble(n))
        enddo
      enddo
      return
      end

