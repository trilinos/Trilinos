      subroutine swaxpy( n, sa, sx, sy, sw )
      implicit none
*
*     constant times a vector plus a vector.
      integer n
      real*8 sa
      real*8 sx( n ), sy( n ), sw( n )
      integer i
c
      do 10 i = 1, n
         sw( i ) = sa*sx( i ) + sy( i )
   10 continue
      return
      end
