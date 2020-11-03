C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDSPLN (STORE, NDB, SLTOP, SLBOT, NSPL,
     *     ZS, XS, YS)
C=======================================================================

      LOGICAL STORE, MATSTR
      REAL    ZS(*), XS(*), YS(*)
      INTEGER NSPL
      REAL    SLTOP(2), SLBOT(2)
      INTEGER NDB

      PARAMETER (BINGO = 1.0E38)

      PARAMETER (MXFLD = 4)
      REAL RVAL(MXFLD)
      INTEGER KVAL(MXFLD), IVAL(MXFLD)
      CHARACTER*8  CVAL(MXFLD)

      REWIND (NDB)
      IPTA = 0
      SLTOP(1) = BINGO
      SLBOT(1) = BINGO
      SLTOP(2) = BINGO
      SLBOT(2) = BINGO

   10 CONTINUE
      CALL FREFLD ( NDB, 0, 'AUTO', MXFLD, IERR,
     *     NFLD, KVAL, CVAL, IVAL, RVAL)
      IF (IERR .EQ. 0) THEN
         IF (KVAL(1) .EQ. 0) THEN
            IF (MATSTR(CVAL(1), 'TOP', 1) .OR.
     &           MATSTR(CVAL(1), 'FRONT', 1)) THEN
            ELSE IF (MATSTR(CVAL(1), 'BOTTOM', 1) .OR.
     &              MATSTR(CVAL(1), 'BACK', 1)) THEN
            ELSE IF (MATSTR(CVAL(1), 'SLOPE', 1)) THEN
               IF (MATSTR(CVAL(2), 'TOP', 1) .OR.
     &              MATSTR(CVAL(2), 'FRONT', 1)) THEN
                  IF (STORE) THEN
                     SLTOP(1) = RVAL(3)
                     SLTOP(2) = RVAL(4)
                  END IF
               ELSE IF (MATSTR(CVAL(2), 'BACK', 1) .OR.
     &                 MATSTR(CVAL(2), 'BOTTOM', 1)) THEN
                  IF (STORE) THEN
                     SLBOT(1) = RVAL(3)
                     SLBOT(2) = RVAL(4)
                  END IF
               END IF
            END IF
         ELSE IF (NFLD .GE. 2) THEN
            IPTA = IPTA + 1
            IF (STORE) THEN
               ZS(IPTA)   = RVAL(1)
               XS(IPTA)   = RVAL(2)
               YS(IPTA)   = RVAL(3)
            END IF
         END IF
         GO TO 10
      END IF

      IF (.NOT. STORE) THEN
         NSPL = IPTA
      END IF

      RETURN
      END
