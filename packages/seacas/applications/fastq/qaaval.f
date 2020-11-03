C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QAAVAL (MXND, NODES, ANGLES, QRAT, AREA, XN, YN, CAREA)
C***********************************************************************

C  SUBROUTINE QAAVAL = CALCULATES THE INTERIOR ANGLES OF A QUAD AND
C                      THE RATIO OF LARGEST TO SMALLEST ANGLE

C***********************************************************************

C  VARIABLES USED:
C     NODES  = THE FOUR NODES OF THE ELEMENT IN CCW ORDER
C     ANGLES = THE FOUR INTERIOR ANGLES IN THE ORDER 4-1, 1-2, 2-3, 3-4
C     QRAT   = THE RATIO OF THE LARGEST TO THE SMALLEST ANGLE
C     CAREA  = .TRUE. IF THE AREA IS TO BE CALCULATED
C     AREA   = ELEMENT'S AREA

C***********************************************************************

      DIMENSION NODES (4), ANGLES (4), AG (4), XN (MXND), YN (MXND)

      LOGICAL CAREA

      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI

      IF (CAREA) THEN
         N1 = NODES (1)
         N2 = NODES (2)
         N3 = NODES (3)
         N4 = NODES (4)
         AREA = 0.5 *  ( (XN (N3) - XN (N1)) *  (YN (N4) - YN (N2))
     &      -  (YN (N3) - YN (N1)) *  (XN (N4) - XN (N2)))
      ENDIF

      DO 100 I = 1, 4
         IF (I .EQ. 4) THEN
            J = 1
         ELSE
            J = I + 1
         ENDIF
         N1 = NODES (I)
         N2 = NODES (J)
         IF ( (XN (N2) .EQ. XN (N1)) .AND. (YN (N2) .EQ. YN (N1))) THEN
            QRAT = 1.E + 10
            RETURN
         ELSE
            AG (J) = ATAN2 (YN (N2) - YN (N1), XN (N2) - XN (N1))
         ENDIF
  100 CONTINUE

      DO 110 J = 1, 4
         IF (J .EQ. 1) THEN
            I = 4
         ELSE
            I = J - 1
         ENDIF
         DIFF = AG (J) - AG (I)
         IF (DIFF .GT. PI) THEN
            DIFF = DIFF - TWOPI
         ELSEIF (DIFF .LT.  - PI) THEN
            DIFF = DIFF + TWOPI
         ENDIF
         ANGLES (J) = PI - DIFF
  110 CONTINUE

      QMIN = ANGLES (1)
      QMAX = ANGLES (1)
      DO 120 I = 2, 4
         QMIN = AMIN1 (QMIN, ANGLES (I))
         QMAX = AMAX1 (QMAX, ANGLES (I))
  120 CONTINUE
      IF (QMIN .GT. 0.) THEN
         QRAT = QMAX / QMIN
      ELSE
         QRAT = 1.0E10
      ENDIF

      RETURN

      END
