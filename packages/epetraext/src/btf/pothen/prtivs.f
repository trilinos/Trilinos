      subroutine   prtivs   ( title, n, x, output )
 
c     ==================================================================
c     ====  bcs sparse matrix package, release 0                    ====
c     ==================================================================
c     ==================================================================
c     ====  prtivs -- print integer vector in table (short format)  ====
c     ==================================================================
c     ==================================================================
c
c     prtivs prints out the integer vector  x  of length  n  to  logical 
c     unit output in a short format.  the character string in  title  is
c     printed as a title for the table
c
c     last modified   --- july 10, 1989   -- jgl --
c
c     --------------
c     ... parameters
c     --------------
 
      character * (*)     title

      integer             n, output
 
      integer             x(n)

c     -------------------
c     ... local variables
c     -------------------

      integer             i, l

      character * 75      line

c     ==================================================================
 
c     ---------------
c     ... write title
c     ---------------
 
      l = min ( len (title), 75 )
 
      do 100 i = 1, l
          line(i:i) = '-'
  100 continue
 
      do 200 i = l+1, 75
          line(i:i) = ' '
  200 continue
 
      write ( output, 2000 ) title (1:l), line (1:l)
 
c     ------------------------
c     ... write out the vector
c     ------------------------
 
      write ( output, 2100 ) x
 
      return
 
c     -----------
c     ... formats
c     -----------

 2000 format ( /5x, a  /5x, a  / )
 
 2100 format ( (10x, 10 (1x, i5)) )

      end
