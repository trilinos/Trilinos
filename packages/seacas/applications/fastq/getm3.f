C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETM3 (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, M1A, M1B,
     &   M2A, M2B, M3A, M3B, XCEN, YCEN, CCW, ERR)
C***********************************************************************

C  SUBROUTINE GETM3 = GETS THE APPROPRIATE SIDE LENGTHS AND DIVISIONS
C                     FOR A TRIANGULAR REGION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS

C***********************************************************************

C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE

C***********************************************************************

      DIMENSION NNPS (MNNPS), ISLIST (NS), LINKL (2, ML), LINKS (MS*2)
      DIMENSION NINT (ML), NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION X (NPER), Y (NPER), NID (NPER), ANGLE (NPER)

      LOGICAL CCW, ERR

C  CALCULATE THE NUMBER OF NODES PER SIDE

      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN
      IF (.NOT. CCW) CALL IREVER (NNPS, NS)

C  FIND THE BEST CORNER NODES IN THE LIST

      CALL PICKM3 (NPER, X, Y, ANGLE, M1, M2, IFIRST)
      IF (IFIRST .NE. 1) CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)

C  NOW SORT THE LIST SO THE LONGEST SIDE IS FIRST

      M3 = NPER - M1 - M2
      MMAX = MAX0 (M1, M2, M3)
      IF (M1 .EQ. MMAX)THEN
         CONTINUE
      ELSEIF (M2 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + 1)
         MHOLD = M1
         M1 = M2
         M2 = M3
         M3 = MHOLD
      ELSEIF (M3 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + M2 + 1)
         MHOLD = M1
         M1 = M3
         M3 = M2
         M2 = MHOLD
      ENDIF

C  SPLIT THE SIDES INTO LOGICAL DIVISIONS

C      IF (M2 .EQ. M3)THEN
C         M1A = (.5 * DBLE(M1)) + .001
C         M1B = M1A
C      ELSEIF (M2 .LT. M3)THEN
C         FACT = DBLE(M2 - 1) / DBLE(M2 + M3 - 2)
C         M1A = FACT * M1 + .001
C         M1A = MAX0 (M1A, 1)
C         M1A = MIN0 (M1A, M1-1)
C         M1B = M1 - M1A
C      ELSE
C         FACT = DBLE(M3 - 1) / DBLE(M2 + M3 - 2)
C         M1B = FACT * M1 + .001
C         M1B = MAX0 (M1B, 1)
C         M1B = MIN0 (M1B, M1-1)
C         M1A = M1-M1B
C      ENDIF
      M1A =  (M1 + M2 - M3) / 2
      M1B = M1 - M1A
      M2A = M2 - M1A
      M2B = M1A
      M3A = M1B
      M3B = M3 - M3A

      ERR = .TRUE.
      IF (M3B .NE. M2A) THEN
         CALL MESAGE ('ERROR GENERATING TRIANGLE DIVISION POINT')
         RETURN
      ENDIF

C  DEFINE THE MIDDLE POINT AS THE AVERAGE OF PROPORIONAL DIVISIONS
C  OF SIDE DIVISION POINT TO OPPOSITE TRIANGLE CORNER LINES

      I1 = 1
      I2 = M1 + 1
      I3 = M1 + M2 + 1
      J1 = I1 + M1A
      J2 = I2 + M2A
      J3 = I3 + M3A

C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS

      D1 = SQRT ( (X (I2) - X (I1)) **2 + (Y (I2) - Y (I1)) **2)
      D2 = SQRT ( (X (I3) - X (I2)) **2 + (Y (I3) - Y (I2)) **2)
      D3 = SQRT ( (X (I1) - X (I3)) **2 + (Y (I1) - Y (I3)) **2)
      D1A = SQRT ( (X (J1) - X (I1)) **2 + (Y (J1) - Y (I1)) **2)
      D1B = SQRT ( (X (I2) - X (J1)) **2 + (Y (I2) - Y (J1)) **2)
      D2A = SQRT ( (X (J2) - X (I2)) **2 + (Y (J2) - Y (I2)) **2)
      D2B = SQRT ( (X (I3) - X (J2)) **2 + (Y (I3) - Y (J2)) **2)
      D3A = SQRT ( (X (J3) - X (I3)) **2 + (Y (J3) - Y (I3)) **2)
      D3B = SQRT ( (X (I1) - X (J3)) **2 + (Y (I1) - Y (J3)) **2)

C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS

      PRO1 = .5 * ( (D3A / D3) + (D1B / D1))
      X1 = X (J2) - (PRO1 * (X (J2) - X (I1)))
      Y1 = Y (J2) - (PRO1 * (Y (J2) - Y (I1)))
      PRO2 = .5 * ( (D2B / D2) + (D1A / D1))
      X2 = X (J3) - (PRO2 * (X (J3) - X (I2)))
      Y2 = Y (J3) - (PRO2 * (Y (J3) - Y (I2)))
      PRO3 = .5 * ( (D2A / D2) + (D3B / D3))
      X3 = X (J1) - (PRO3 * (X (J1) - X (I3)))
      Y3 = Y (J1) - (PRO3 * (Y (J1) - Y (I3)))

C  AVERAGE POINTS TO GET THE CENTER

      XCEN =  (X1 + X2 + X3) / 3.
      YCEN =  (Y1 + Y2 + Y3) / 3.

      ERR = .FALSE.
      RETURN

      END
