C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PLANE3 (COORD, NUMNP, DIST, DISTR, NDIM, P1, P2, TOLER,
     *   NODEL, SORTYP, MAP, SORUP, INUM, OPT, SELECT)
      DIMENSION COORD (NUMNP,*), DIST(*), DISTR(*), P1(*), P2(*),
     *   TOLER(*), MAP(*)
      CHARACTER*(*) NODEL, SORTYP, OPT
      LOGICAL SORUP, SELECT(*), ISABRT
      include 'nu_io.blk'

      CALL LOCOUT ('PLANE', NDIM, NODEL, TOLER, SORTYP, P1, P2, ' ')

      A = P2(1)
      B = P2(2)
      C = P2(3)
      D = A * P1(1) + B * P1(2) + C * P1(3)

      TEMP = TOLER(1)
      TOLER(1) = MAX(0.0, TEMP - TOLER(2))
      TOLER(2) = MAX(0.0, TEMP + TOLER(2))

      DO 10 I=1, NUMNP
         IF (SELECT(I)) THEN
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            Z0 = COORD(I,3)
            DIST(I) = ABS(A * X0 + B * Y0 + C * Z0 - D) /
     *         SQRT(A**2 + B**2 + C**2)
            DISTR(I)=(X0 - P1(1))**2 + (Y0 - P1(2))**2 + (Z0 - P1(3))**2
         END IF
   10 CONTINUE

      INUM = 0
      DISMIN = 1.0E30
      DO 20 I=1, NUMNP
         IF (SELECT(I)) THEN
            DISMIN = MIN(DIST(I), ABS(DISMIN-TEMP))
            IF (DIST(I) .GE. TOLER(1) .AND. DIST(I) .LE. TOLER(2)) THEN
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
      ELSE IF (SORTYP .EQ. 'RADIAL') THEN
         CALL INDEXX (DISTR,      MAP, INUM, .FALSE.)
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
   30    CONTINUE
   40    FORMAT (/,T50,'DISTANCE',/2X,A8,T16,'X',T26,'Y',T36,'Z',
     *      T45,'NORMAL',T57,'RADIAL',/)
         DO 60 IN = IBEG, IEND, IINC
            IF (ISABRT()) RETURN
            I = MAP(IN)
            DO 50 IO=IOMIN, IOMAX
               WRITE (IO, 90) I, (COORD(I,J),J=1,3), DIST(I),
     *            SQRT(DISTR(I))
   50       CONTINUE
   60    CONTINUE
         IF (INUM .EQ. 0) THEN
            DO 70 IO=IOMIN, IOMAX
               WRITE (IO, 80) DISMIN
   70       CONTINUE
         END IF
      END IF
   80 FORMAT (/' None found within tolerance, minimum distance = ',
     *   1PE12.3,/)
   90 FORMAT (I10, 3(F10.4), 2(1PE12.3))
      RETURN
      END
