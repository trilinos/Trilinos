C $Id: getm3.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: getm3.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1998/07/14 18:19:01  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:08:26  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:08:24  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]GETM3.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETM3 (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, M1A, M1B,
     &   M2A, M2B, M3A, M3B, XCEN, YCEN, CCW, ERR)
C***********************************************************************
C
C  SUBROUTINE GETM3 = GETS THE APPROPRIATE SIDE LENGTHS AND DIVISIONS
C                     FOR A TRIANGULAR REGION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE
C
C***********************************************************************
C
      DIMENSION NNPS (MNNPS), ISLIST (NS), LINKL (2, ML), LINKS (MS*2)
      DIMENSION NINT (ML), NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION X (NPER), Y (NPER), NID (NPER), ANGLE (NPER)
C
      LOGICAL CCW, ERR
C
C  CALCULATE THE NUMBER OF NODES PER SIDE
C
      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN
      IF (.NOT. CCW) CALL IREVER (NNPS, NS)
C
C  FIND THE BEST CORNER NODES IN THE LIST
C
      CALL PICKM3 (NPER, X, Y, ANGLE, M1, M2, IFIRST)
      IF (IFIRST .NE. 1) CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)
C
C  NOW SORT THE LIST SO THE LONGEST SIDE IS FIRST
C
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
C
C  SPLIT THE SIDES INTO LOGICAL DIVISIONS
C
C      IF (M2 .EQ. M3)THEN
C         M1A = (.5 * FLOAT (M1)) + .001
C         M1B = M1A
C      ELSEIF (M2 .LT. M3)THEN
C         FACT = FLOAT (M2 - 1) / FLOAT (M2 + M3 - 2)
C         M1A = FACT * M1 + .001
C         M1A = MAX0 (M1A, 1)
C         M1A = MIN0 (M1A, M1-1)
C         M1B = M1 - M1A
C      ELSE
C         FACT = FLOAT (M3 - 1) / FLOAT (M2 + M3 - 2)
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
C
      ERR = .TRUE.
      IF (M3B .NE. M2A) THEN
         CALL MESAGE ('ERROR GENERATING TRIANGLE DIVISION POINT')
         RETURN
      ENDIF
C
C  DEFINE THE MIDDLE POINT AS THE AVERAGE OF PROPORIONAL DIVISIONS
C  OF SIDE DIVISION POINT TO OPPOSITE TRIANGLE CORNER LINES
C
      I1 = 1
      I2 = M1 + 1
      I3 = M1 + M2 + 1
      J1 = I1 + M1A
      J2 = I2 + M2A
      J3 = I3 + M3A
C
C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS
C
      D1 = SQRT ( (X (I2) - X (I1)) **2 + (Y (I2) - Y (I1)) **2)
      D2 = SQRT ( (X (I3) - X (I2)) **2 + (Y (I3) - Y (I2)) **2)
      D3 = SQRT ( (X (I1) - X (I3)) **2 + (Y (I1) - Y (I3)) **2)
      D1A = SQRT ( (X (J1) - X (I1)) **2 + (Y (J1) - Y (I1)) **2)
      D1B = SQRT ( (X (I2) - X (J1)) **2 + (Y (I2) - Y (J1)) **2)
      D2A = SQRT ( (X (J2) - X (I2)) **2 + (Y (J2) - Y (I2)) **2)
      D2B = SQRT ( (X (I3) - X (J2)) **2 + (Y (I3) - Y (J2)) **2)
      D3A = SQRT ( (X (J3) - X (I3)) **2 + (Y (J3) - Y (I3)) **2)
      D3B = SQRT ( (X (I1) - X (J3)) **2 + (Y (I1) - Y (J3)) **2)
C
C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS
C
      PRO1 = .5 * ( (D3A / D3) + (D1B / D1))
      X1 = X (J2) - (PRO1 * (X (J2) - X (I1)))
      Y1 = Y (J2) - (PRO1 * (Y (J2) - Y (I1)))
      PRO2 = .5 * ( (D2B / D2) + (D1A / D1))
      X2 = X (J3) - (PRO2 * (X (J3) - X (I2)))
      Y2 = Y (J3) - (PRO2 * (Y (J3) - Y (I2)))
      PRO3 = .5 * ( (D2A / D2) + (D3B / D3))
      X3 = X (J1) - (PRO3 * (X (J1) - X (I3)))
      Y3 = Y (J1) - (PRO3 * (Y (J1) - Y (I3)))
C
C  AVERAGE POINTS TO GET THE CENTER
C
      XCEN =  (X1 + X2 + X3) / 3.
      YCEN =  (Y1 + Y2 + Y3) / 3.
C
      ERR = .FALSE.
      RETURN
C
      END
