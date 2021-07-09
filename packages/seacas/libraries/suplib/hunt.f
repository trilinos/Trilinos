C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      subroutine hunt(a,n,x,jlo)

C ... find jlo such that a(jlo) .le. x .and. a(jlo+1) .gt. x
C     (if jlo .ne. n)

C     Start search at passed in 'jlo' position

      DIMENSION a(n)

      integer jlo, low, high

      if (jlo .lt. 1 .or. jlo .gt. n) jlo = 1

      if (a(jlo) .eq. x) then
        return
      else if (a(1) .ge. x) then
        jlo = 1
        return
      else if (a(n) .lt. x) then
        jlo = n
        return
      end if

      if (a(jlo) .le. x) then
        low  = jlo
        high = n-1
      else
        low  = 1
        high = jlo
      end if

      do jlo = low, high
        if (a(jlo) .le. x .and. a(jlo+1) .ge. x) then
          return
        end if
      end do
C ... should not get here since extremes checked above...

      end
