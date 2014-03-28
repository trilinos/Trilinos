C $Id: vol2d.f,v 1.2 1992/12/11 22:34:15 gdsjaar Exp $
C $Log: vol2d.f,v $
C Revision 1.2  1992/12/11 22:34:15  gdsjaar
C Fixed problem with incorrect determination of cavity center in 2d
C
c Revision 1.1.1.1  1991/02/21  15:46:17  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:46:16  gdsjaar
c Initial revision
c
      SUBROUTINE VOL2D( COORD, LSTSN, NSEG, VOLUME, NDIM, NUMESS, AXI,
     *   CENT, NUMNP, CENTER)
C
C***********************************************************************
C
C     DESCRIPTION:
C       This routine computes the volume of a cavity formed
C       by the boundary of an element side set flag
C
C     FORMAL PARAMETERS:
C       COORD   REAL      Nodal Coordinates
C       LSTSN   INTEGER   List of nodes on this boundary
C       NSEG    INTEGER   Number of segments in this boundary
C       VOLUME  REAL      Volume of this cavity
C       NDIM    INTEGER   Number of Dimensions
C       NUMESS  INTEGER   Number of Element Side Set Flags
C       AXI     LOGICAL   TRUE if axisymmetric mesh
C       CENT    REAL      Apex of volume triangles
C
C     CALLED BY:
C
C***********************************************************************
C
      LOGICAL AXI, CENTER
      DIMENSION COORD(NUMNP, *), LSTSN(*), CENT(2)
      TWOPI = 2.0 * ATAN2(0.0, -1.0)
C
      VOLUME = 0.0
C
C ... CALCULATE APPROXIMATE CENTROID OF CAVITY.  ASSUME XC ON
C     SYMMETRY AXIS.
C
      IF (.NOT. CENTER) THEN
         DO 10 KSEG = 1 , 2*NSEG
             J = LSTSN(KSEG)
C        
             CENT(2) = CENT(2) + COORD(J,2)
   10    CONTINUE
         CENT(2) = CENT(2) / (2 * NSEG)
      END IF
C
      DO 20 KSEG = 1 , NSEG
          J = LSTSN(2*KSEG)
          I = LSTSN(2*KSEG - 1)
C
          X1 = COORD(I,1) - CENT(1)
          X2 = COORD(J,1) - CENT(1)
C
          Y1 = COORD(I,2) - CENT(2)
          Y2 = COORD(J,2) - CENT(2)
C
          VP = (Y1 * X2 - Y2 * X1) / 2.0
          IF (AXI) THEN
              XC = (X1 + X2) / 3.0
              VP = TWOPI * XC * VP
          END IF
          VOLUME = VOLUME + VP
   20 CONTINUE
C
      RETURN
      END
