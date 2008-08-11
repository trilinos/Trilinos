      subroutine   izero   ( n, x, incx )
 
c     ==================================================================
c     ==================================================================
c     ====  izero -- initialize integer vector to zero              ====
c     ==================================================================
c     ==================================================================
 
c     purpose ... initializes integer vector to zero
 
c     created       ... mar. 8, 1985
c     last modified ... apr. 19, 1985
 
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             n, incx
 
      integer             x (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             xaddr, i
 
c     ==================================================================
 
      if  ( incx .eq. 1 )  then
 
c         ----------------------------------
c         ... unit increment (standard case)
c         ----------------------------------
 
          do 100 i = 1, n
              x(i) = 0
  100     continue
 
      else
 
c         ----------------------
c         ... non-unit increment
c         ----------------------
 
          xaddr = 1
          if  ( incx .lt. 0 )  then
              xaddr = (-n+1)*incx + 1
          endif
 
          do 200 i = 1, n
              x (xaddr) = 0
              xaddr     = xaddr + incx
  200     continue
 
      endif
 
      return
 
      end
