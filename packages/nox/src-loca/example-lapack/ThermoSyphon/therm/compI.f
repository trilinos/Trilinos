c
c subroutine compI to compute matrix Identity 
c
      subroutine compI(m,matI, mmax)
      implicit none 
      integer m, i, j, mmax
      double precision matI(0:mmax, 0:mmax)
      do i=0, mmax
        do j=0, mmax
          matI(i,j)=0.d0
        enddo
      enddo
c
      do i=0,m
       matI(i,i)=1.d0
      enddo
      return
      end

