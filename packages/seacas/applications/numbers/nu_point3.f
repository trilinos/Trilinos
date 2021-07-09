C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE POINT3 (COORD, NUMNP, DIST, NDIM, P1, TOLER,
     *   NODEL, SORTYP, MAP, ANGLE, SORUP, INUM, OPT, SELECT)
      DIMENSION COORD (NUMNP,*), DIST(*), P1(*), TOLER(*),
     *   MAP(*), ANGLE(2,*)
      CHARACTER*(*) NODEL, SORTYP, OPT
      LOGICAL SORUP, SELECT(*), ISABRT
      include 'nu_io.blk'
      PI = ATAN2(0.0, -1.0)

      CALL LOCOUT ('POINT', NDIM, NODEL, TOLER, SORTYP, P1, P1, ' ')

      TEMP = TOLER(1)
      TOLER(1) = MAX(0.0, TEMP - TOLER(2))
      TOLER(2) = MAX(0.0, TEMP + TOLER(2))

      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)

      DO 10 I=1, NUMNP
         IF (SELECT(I)) THEN
         X0 = COORD(I,1)
         Y0 = COORD(I,2)
         Z0 = COORD(I,3)

         DIST(I) = (X1 - X0)**2 + (Y1 - Y0)**2 + (Z1 - Z0)**2
         END IF
   10 CONTINUE

      INUM = 0
      DISMIN = 1.0E38
      DO 20 I=1, NUMNP
         IF (SELECT(I)) THEN
         DISMIN = MIN(DIST(I), ABS(DISMIN-TEMP))
         IF (DIST(I) .GE. TOLER(1)**2 .AND. DIST(I) .LE. TOLER(2)**2)
     *      THEN
            DX = COORD(I,1) - P1(1)
            DY = COORD(I,2) - P1(2)
            DZ = COORD(I,3) - P1(3)
            DS = SQRT(DX**2 + DZ**2)

            FIX1 = SIGN(0.5,ABS(DX+DZ)) + SIGN(0.5,-ABS(DX+DZ))
            FIX2 = SIGN(0.5,ABS(DS+DY)) + SIGN(0.5,-ABS(DS+DY))
            ANGLE(1,I) = ATAN2(DZ,DX+FIX1) * 180.0 / PI
            ANGLE(2,I) = ATAN2(DY,DS+FIX2) * 180.0 / PI

            INUM = INUM + 1
            MAP(INUM) = I
         END IF
         END IF
   20 CONTINUE

      IF      (SORTYP .EQ. 'X') THEN
         CALL INDEXX (COORD(1,1), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Y') THEN
         CALL INDEXX (COORD(1,2), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Z') THEN
         CALL INDEXX (COORD(1,3), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'THETA') THEN
         CALL INDEXN (ANGLE,      2, 1, MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'PHI') THEN
         CALL INDEXN (ANGLE,      2, 2, MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'DISTANCE') THEN
         CALL INDEXX (DIST,       MAP, INUM, .FALSE.)
      END IF

      IF (SORUP) THEN
         IBEG = 1
         IEND = INUM
         IINC = 1
      ELSE
         IBEG = INUM
         IEND = 1
         IINC = -1
      END IF

      IF (OPT .EQ. '*' .OR. INDEX(OPT, 'P') .GT. 0) THEN
      DO 30 IO=IOMIN, IOMAX
         WRITE (IO, 40) NODEL
   30 CONTINUE
   40 FORMAT (/,2X,A8,T16,'X',T26,'Y',T36,'Z',T44,'DISTANCE',
     *   T58,'THETA',T71,'PHI')
      DO 60 IN = IBEG, IEND, IINC
         IF (ISABRT()) RETURN
         I = MAP(IN)
         DO 50 IO=IOMIN, IOMAX
            WRITE (IO, 90) I, (COORD(I,J),J=1,3), SQRT(DIST(I)),
     *         ANGLE(1,I), ANGLE(2,I)
   50    CONTINUE
   60 CONTINUE

      IF (INUM .EQ. 0) THEN
         DO 70 IO=IOMIN, IOMAX
            WRITE (IO, 80) SQRT(DISMIN)
   70    CONTINUE
      END IF
      END IF
   80 FORMAT (/' None found within range, minimum distance = ',
     *   1PE12.3,/)
   90 FORMAT (I10, 3(F10.4), 3(1PE12.3))
      RETURN
      END
