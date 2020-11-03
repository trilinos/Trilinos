C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETM1(ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, SCHSTR,
     &   M1, CCW, NORM, REAL, ERR)
C***********************************************************************

C  SUBROUTINE GETM1 = GETS THE APPROPRIATE M1 VALUE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS

C***********************************************************************

C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE

C***********************************************************************

      DIMENSION NNPS(MNNPS), ISLIST(NS), LINKL(2, ML), LINKS(MS*2)
      DIMENSION NINT(ML), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION X(NPER), Y(NPER), NID(NPER), ANGLE(NPER)

      LOGICAL CCW, ERR, NORM, REAL

      CHARACTER*72 SCHSTR

C  SEE IF LETTER M OCCURS IN THE SCHEME (NOT THE "NORMAL" CASE)

      NORM = .TRUE.
      DO 100 J = 1, 72
         IF ( (SCHSTR (J:J) .EQ. 'M') .OR.
     &      (SCHSTR (J:J) .EQ. 'm') ) THEN
            NORM = .FALSE.
            GOTO 110
         ENDIF
  100 CONTINUE
  110 CONTINUE

C  CALCULATE THE NUMBER OF NODES PER SIDE

      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN

C  NORMAL CASE - TRY THE FIRST SIDE IN THE SIDE LIST AS THE M1 BASE
C  IN CCW ORDER IF IT IS NOT RIDICULOUS

      IF (NORM) THEN
         M1 = NNPS (1) - 1
         IF (.NOT. CCW) M1 = NPER / 2 - M1
         IF ( (M1 .GT. 0).AND. (M1 .LE. NPER/2) ) RETURN
      ENDIF
      IF (.NOT. CCW) CALL IREVER (NNPS, NS)

C  IF THE BOUNDARY IS A LOGICAL RECTANGLE,  USE IT

      IF ( (NS .EQ. 4) .AND. (NNPS (1) .EQ. NNPS (3)) .AND.
     &   (NNPS (2) .EQ. NNPS (4) ) ) THEN
         M1 = NNPS (1) - 1

C  OTHERWISE,  FIND AN INITIAL M1 FOR A NON-LOGICAL RECTANGLE

      ELSE
         CALL PICKM1 (NPER, X, Y, ANGLE, M1, IFIRST, REAL)
         IF (IFIRST .NE. 1) CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)
      ENDIF

      RETURN

      END
