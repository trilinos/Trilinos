
cSubroutine compR to compute matrix R
      subroutine compR(m,matR, mmax)
      implicit none 
      integer m, mmax, i, j
      double precision matR(0:mmax,0:mmax)
      do i=0, mmax
        do j=0, mmax
          matR(i,j)=0.d0
        enddo
       enddo
       do i=2,m
         matR(i,i-1)=.5d0
         matR(i-1,i)=.5d0
       enddo 
      matR(1,0)=1.d0
      matR(1,2)=.5d0
      return
      end

