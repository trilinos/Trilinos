C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RANGE (LEN, LIST, IOMIN, IOMAX)
      LOGICAL LIST(*), INRNG
      INTEGER IRANGE(3)
      CHARACTER*80 LINE

      INRNG = .FALSE.
      IRINC = 0

      DO 10 I=1, LEN
         IF (LIST(I)) THEN
            IBEG = I+1
            LASTSL = I
            GO TO 20
         END IF
   10 CONTINUE

   20 CONTINUE
      DO 40 I = IBEG, LEN
         IF (LIST(I)) THEN
            IRINC1 = I - LASTSL
            IF (IRINC .EQ. 0) IRINC = IRINC1
            IF (.NOT. INRNG) THEN
               IRBEG = LASTSL
               IRINC = IRINC1
               INRNG = .TRUE.
            ELSE IF (INRNG .AND. IRINC1 .NE. IRINC) THEN
               IREND = I
               INRNG = .FALSE.
               IRANGE(1) = IRBEG
               IRANGE(2) = IREND - IRINC1
               IRANGE(3) = IRINC
               LINE = ' '
               CALL FFADDV (IRANGE, LINE)
               DO 30 IO = IOMIN, IOMAX
                  CALL LOGERR ('CMDSPEC', LINE(:LENSTR(LINE)), IO)
   30          CONTINUE
               IRINC = IRINC1
            END IF
            LASTSL = I
         END IF
   40 CONTINUE
      IRANGE(1) = IRBEG
      IRANGE(2) = LASTSL
      IRANGE(3) = IRINC
      LINE = ' '
      CALL FFADDV (IRANGE, LINE)
      DO 50 IO = IOMIN, IOMAX
         CALL LOGERR ('CMDSPEC', LINE(:LENSTR(LINE)), IO)
   50 CONTINUE
      RETURN
      END
