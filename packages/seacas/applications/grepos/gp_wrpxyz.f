C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRPXYZ (X, Y, Z, NUMNP, IWARP, NORMAL, REFDIS)
C=======================================================================
      REAL X(NUMNP), Y(NUMNP), Z(NUMNP)

C ... Origin
      if (iwarp .eq. 1) then
        call PRTERR('PROGRAM', 'Origin warping not implemented yet.')
C ... Xaxis
      else if (iwarp .eq. -1) then
        if (normal .eq. 2) then
          call warpit(x, y, z, numnp, refdis)
        else if (normal .eq. 3) then
          call warpit(x, z, y, numnp, refdis)
        end if

C ... Yaxis
      else if (iwarp .eq. -2) then
        if (normal .eq. 3) then
          call warpit(y, z, x, numnp, refdis)
        else if (normal .eq. 1) then
          call warpit(y, x, z, numnp, refdis)
        end if

C ... Zaxis
      else if (iwarp .eq. -3) then
        if (normal .eq. 1) then
          call warpit(z, x, y, numnp, refdis)
        else if (normal .eq. 2) then
          call warpit(z, y, x, numnp, refdis)
        end if
      end if

      RETURN
      END

      SUBROUTINE WARPIT(C1, C2, C3, NUMNP, REFDIS)
      REAL C1(NUMNP), C2(NUMNP), C3(NUMNP)

      do 10 i=1, numnp
        c1(i) = c1(i)

        radius = c2(i)
        theta  = c3(i) / refdis

        c3(i) = radius * sin(theta)
        c2(i) = radius * cos(theta)
 10   continue

      return
      end
