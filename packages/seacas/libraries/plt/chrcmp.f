C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION CHRCMP(KWD,PART1,PART2)
      CHARACTER*(*) KWD
      CHARACTER*(*) PART1
      CHARACTER*(*) PART2

      CALL CHRTRM(KWD,LK)
      CALL CHRTRM(PART1,LF)
      CALL CHRTRM(PART2,LV)
      IF (LK.LT.LF .OR. LK.GT.LF+LV) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (KWD(1:LF).NE.PART1(1:LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (LK.EQ.LF) THEN
         CHRCMP = .TRUE.
         RETURN

      END IF

      IF (KWD(LF+1:LK).NE.PART2(1:LK-LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      CHRCMP = .TRUE.
      RETURN

      END
