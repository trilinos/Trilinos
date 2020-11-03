C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE VOL2D( COORD, LSTSN, NSEG, VOLUME, NDIM, NUMESS, AXI,
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
C       NDIM    INTEGER   Number of Dimensions
C       NUMESS  INTEGER   Number of Element Side Set Flags
C       AXI     LOGICAL   TRUE if axisymmetric mesh
C       CENT    REAL      Apex of volume triangles

C     CALLED BY:

C***********************************************************************

      LOGICAL AXI, CENTER
      DIMENSION COORD(NUMNP, *), LSTSN(*), CENT(2)
      TWOPI = 2.0 * ATAN2(0.0, -1.0)

      VOLUME = 0.0

C ... CALCULATE APPROXIMATE CENTROID OF CAVITY.  ASSUME XC ON
C     SYMMETRY AXIS.

      IF (.NOT. CENTER) THEN
         DO 10 KSEG = 1 , 2*NSEG
             J = LSTSN(KSEG)

             CENT(2) = CENT(2) + COORD(J,2)
   10    CONTINUE
         CENT(2) = CENT(2) / (2 * NSEG)
      END IF

      DO 20 KSEG = 1 , NSEG
          J = LSTSN(2*KSEG)
          I = LSTSN(2*KSEG - 1)

          X1 = COORD(I,1) - CENT(1)
          X2 = COORD(J,1) - CENT(1)

          Y1 = COORD(I,2) - CENT(2)
          Y2 = COORD(J,2) - CENT(2)

          VP = (Y1 * X2 - Y2 * X1) / 2.0
          IF (AXI) THEN
              XC = (X1 + X2) / 3.0
              VP = TWOPI * XC * VP
          END IF
          VOLUME = VOLUME + VP
   20 CONTINUE

      RETURN
      END
