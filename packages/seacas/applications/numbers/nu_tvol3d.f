C    Copyright(C) 2025 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TVOL3D( COORD, LSTSN, NSEG, VOLUME, NDIM, NUMESS,
     *   CENT, NUMNP, CENTER)

C***********************************************************************

C     DESCRIPTION:
C       This routine computes the volume of a cavity formed
C       by the boundary of an element side set flag

C     FORMAL PARAMETERS:
C       COORD   REAL      Nodal Coordinates
C       LSTSN   INTEGER   List of nodes on this boundary
C       NSEG    INTEGER   Number of segments in this boundary
C       VOLUME  REAL      Volume of this cavity
C       NDIM    INTEGER   Number of Nodes
C       CENT    REAL      Apex of Cavity Volume pentahedra

C     CALLED BY:

C***********************************************************************

      DIMENSION COORD(NUMNP, *), LSTSN(*), CENT(3)
      LOGICAL CENTER
      REAL M(3,3)
      
      VOLUME = 0.0

      IF (.NOT. CENTER) THEN
        XC = 0.0
        YC = 0.0
        ZC = 0.0
        DO KSEG = 1 , NSEG
          K = LSTSN(3*KSEG)
          J = LSTSN(3*KSEG - 1)
          I = LSTSN(3*KSEG - 2)

          X1 = COORD(I,1)
          X2 = COORD(J,1)
          X3 = COORD(K,1)

          Y1 = COORD(I,2)
          Y2 = COORD(J,2)
          Y3 = COORD(K,2)

          Z1 = COORD(I,3)
          Z2 = COORD(J,3)
          Z3 = COORD(K,3)

          XC = XC + x1 + x2 + x3
          YC = YC + y1 + y2 + y3
          ZC = ZC + z1 + z2 + z3
        END DO
        CENT(1) = XC / (3*NSEG)
        CENT(2) = YC / (3*NSEG)
        CENT(3) = ZC / (3*NSEG)
      END IF

      X4 = CENT(1)
      Y4 = CENT(2)
      Z4 = CENT(3)

      DO 100 KSEG = 1 , NSEG
         K = LSTSN(3*KSEG)
         J = LSTSN(3*KSEG - 1)
         I = LSTSN(3*KSEG - 2)

         X1 = COORD(I,1)
         X2 = COORD(J,1)
         X3 = COORD(K,1)

         Y1 = COORD(I,2)
         Y2 = COORD(J,2)
         Y3 = COORD(K,2)

         Z1 = COORD(I,3)
         Z2 = COORD(J,3)
         Z3 = COORD(K,3)

         m(1,1) = x1 - x2
         m(1,2) = y1 - y2
         m(1,3) = z1 - z2

         m(2,1) = x2 - x3
         m(2,2) = y2 - y3
         m(2,3) = z2 - z3

         m(3,1) = x3 - x4
         m(3,2) = y3 - y4
         m(3,3) = z3 - z4

         det = (m(1,1) * (m(2,2) * m(3,3) - m(2,3) * m(3,2)) -
     $          m(2,1) * (m(1,2) * m(3,3) - m(1,3) * m(3,2)) +
     $          m(3,1) * (m(1,2) * m(2,3) - m(1,3) * m(2,2)))

         VOLUME = VOLUME - det / 6.0
  100 CONTINUE

      RETURN
      END
