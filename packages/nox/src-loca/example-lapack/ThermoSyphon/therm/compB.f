c
c subroutine compB to compute matrix B
c
      subroutine compB(m,matB, mmax)
      implicit none 
      integer m, i, j, mmax
      double precision matB(0:mmax, 0:mmax)
      do i=0, mmax
        do j=0, mmax
          matB(i,j)=0.d0
        enddo
      enddo
c
c  diag zero
c
      do i=2,m
        matB(i,i-1) = 5.d-1/i
        matB(i-1,i) = -5.d-1/(i-1)
      enddo
      matB(1,0)=1.d0
      matB(1,2)=-.5d0
      return
      end

