      double precision function norm ( n, sx )
      implicit none
*
*     Euclidean norm of a real vector
      integer n
      real*8 sx( n )
      integer i
c
      norm = 0.0D0
c
      do 10 i = 1, n
        norm = norm + sx( i )*sx( i )
 10   continue
c
      norm = sqrt(norm)
c
      return
      end
