C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETM5 (IA, ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE,
     &   NLPS, ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER,
     &   M1A, M1B, M2, M3A, M3B, M4A, M4B, M5, MC, XCEN, YCEN, CCW, ERR)
C***********************************************************************

C  SUBROUTINE GETM5 = GETS THE APPROPRIATE SIDE LENGTHS AND DIVISIONS
C                      FOR A PENTAGON REGION

C  WRITTEN BY: HORACIO RECALDE                   DATE: JAN 1988

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS

C***********************************************************************

C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE

C***********************************************************************

      DIMENSION IA(1)
      DIMENSION NNPS(MNNPS), ISLIST(NS), LINKL(2, ML), LINKS(MS*2)
      DIMENSION NLPS(MS), NINT(ML), IFLINE(MS), ILLIST(MS*3)
      DIMENSION X(NPER), Y(NPER), NID(NPER), ANGLE(NPER)
      DIMENSION XJ(3), YJ(3)

      LOGICAL CCW, ERR

C  CALCULATE THE NUMBER OF NODES PER SIDE

      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR)RETURN
      IF (.NOT.CCW) CALL IREVER (NNPS, NS)

C   RESERVE MEMORY FOR THE STACKS

      CALL MDRSRV ('IST2', IP2, NPER)
      CALL MDRSRV ('IST3', IP3, NPER)
      CALL MDRSRV ('IST4', IP4, NPER)
      CALL MDRSRV ('IST5', IP5, NPER)
      CALL MDRSRV ('INDST', INDP, NPER)

C  FIND THE BEST CORNER NODES IN THE LIST

      CALL PICKM5 (NPER, X, Y, ANGLE, IA(IP2), IA(IP3), IA(IP4),
     &   IA(IP5), IA(INDP), IFIRST, M1, M2, M3, M4)
      IF (IFIRST .EQ. 0) THEN
         ERR = .TRUE.
         CALL MESAGE ('ERROR FITTING LOGICAL PENTAGON TO DATA')
         RETURN
      ELSE IF (IFIRST .EQ. -1) THEN
         ERR = .TRUE.
         CALL MESAGE ('TOLERANCE EXCEEDED')
         RETURN
      ELSE IF (IFIRST.NE.1) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)
      END IF

C  DELETE THE STACKS

      CALL MDDEL ('IST2')
      CALL MDDEL ('IST3')
      CALL MDDEL ('IST4')
      CALL MDDEL ('IST5')
      CALL MDDEL ('INDST')

C  NOW SORT THE LIST SO THE LONGEST SIDE IS FIRST

      M5 = NPER - M1 - M2 - M3 - M4
      MMAX = MAX0(M1, M2, M3, M4, M5)
      IF (M1 .EQ. MMAX) THEN
         KNUM = 0
      ELSE IF (M2 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + 1)
         KNUM = 1
      ELSE IF (M3 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + M2 + 1)
         KNUM = 2
      ELSE IF (M4 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + M2 + M3 + 1)
         KNUM = 3
      ELSE IF (M5 .EQ. MMAX) THEN
         CALL FQ_ROTATE (NPER, X, Y, NID, M1 + M2 + M3 + M4 + 1)
         KNUM = 4
      END IF
      DO 100 KK = 1, KNUM
         MHOLD = M1
         M1 = M2
         M2 = M3
         M3 = M4
         M4 = M5
         M5 = MHOLD
  100 CONTINUE

C  SPLIT THE SIDES INTO LOGICAL DIVISIONS

      M1A = (M1 + M4 + M5 - M2 - M3)/2
      M1B = (M1 + M2 + M3 - M4 - M5)/2
      M3A = M1B
      M3B = (M3 + M4 + M5 - M1 - M2)/2
      M4A = (M2 + M3 + M4 - M1 - M5)/2
      M4B = M1A
      MC = (M1 + M2 + M5 - M3 - M4)/2

C  DEFINE THE MIDDLE POINT AS THE AVERAGE OF PROPORIONAL DIVISIONS
C  OF SIDE DIVISION POINT TO OPPOSITE TRIANGLE CORNER LINES

      I1 = M1A + 1
      I2 = I1 + M1B + M2 + M3A
      I3 = I2 + M3B + M4A

C  FIND DISTANCES FROM CORNER TO CORNER, AND CORNERS TO SIDE DIVISIONS

      D1 = SQRT((X(I2) - X(I1))**2 + (Y(I2) - Y(I1))**2)
      D2 = SQRT((X(I3) - X(I2))**2 + (Y(I3) - Y(I2))**2)
      D3 = SQRT((X(I1) - X(I3))**2 + (Y(I1) - Y(I3))**2)
      D1A = DBLE(M4B)*D1/DBLE(M4)
      D1B = D1 - D1A
      D2A = DBLE(M1B)*D2/DBLE(M1)
      D2B = D2 - D1A
      D3A = DBLE(M3B)*D3/DBLE(M3)
      D3B = D3 - D3A
      XJ(1) = X(I1) + (X(I2) - X(I1))*D1A/D1
      YJ(1) = Y(I1) + (Y(I2) - Y(I1))*D1A/D1
      XJ(2) = X(I2) + (X(I3) - X(I2))*D2A/D2
      YJ(2) = Y(I2) + (Y(I3) - Y(I2))*D2A/D2
      XJ(3) = X(I3) + (X(I1) - X(I3))*D3A/D3
      YJ(3) = Y(I3) + (Y(I1) - Y(I3))*D3A/D3

C  GET MIDPOINT TRIALS 1, 2, AND 3 AS PROPORTIONS

      PRO1 = .5*((D3A/D3) + (D1B/D1))
      X1 = XJ(2) - (PRO1*(XJ(2) - X(I1)))
      Y1 = YJ(2) - (PRO1*(YJ(2) - Y(I1)))
      PRO2 = .5*((D2B/D2) + (D1A/D1))
      X2 = XJ(3) - (PRO2*(XJ(3) - X(I2)))
      Y2 = YJ(3) - (PRO2*(YJ(3) - Y(I2)))
      PRO3 = .5*((D2A/D2) + (D3B/D3))
      X3 = XJ(1) - (PRO3*(XJ(1) - X(I3)))
      Y3 = YJ(1) - (PRO3*(YJ(1) - Y(I3)))

C  AVERAGE POINTS TO GET THE CENTER

      XCEN = (X1 + X2 + X3)/3.
      YCEN = (Y1 + Y2 + Y3)/3.

      ERR = .FALSE.
      RETURN

      END
