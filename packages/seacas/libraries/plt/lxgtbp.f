C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXGTBP(STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) STR
      CHARACTER CH
      LOGICAL QFLAG
      CHARACTER QCH

      NS = 0
      CH = ILINE(JLINE:JLINE)
      LXGTBP = (CH.EQ.'(') .OR. (CH.EQ.'[') .OR. (CH.EQ.'{')
      IF (.NOT.LXGTBP) THEN
         RETURN

      END IF

      J = JLINE
      LP = 1
      QFLAG = .FALSE.
 2440 CONTINUE
      J = J + 1
      CH = ILINE(J:J)
      IF (.NOT.QFLAG) THEN
         IF (CH.EQ.'(' .OR. CH.EQ.'[' .OR. CH.EQ.'{') THEN
            LP = LP + 1

         ELSE IF (CH.EQ.')' .OR. CH.EQ.']' .OR. CH.EQ.'}') THEN
            LP = LP - 1

         ELSE IF (CH.EQ.'''' .OR. CH.EQ.'"' .OR. CH.EQ.CHAR(96)) THEN
            QFLAG = .TRUE.
            QCH = CH

         ELSE IF (CH.EQ.CHAR(0)) THEN
            LXGTBP = .FALSE.
            RETURN

         END IF

         IF (LP.EQ.0) THEN
            GO TO 2460

         END IF

         NS = NS + 1
         STR(NS:NS) = CH

      ELSE
         NS = NS + 1
         STR(NS:NS) = CH
         IF (CH.EQ.QCH) THEN
            IF (CH.EQ.ILINE(J+1:J+1)) THEN
               NS = NS + 1
               STR(NS:NS) = CH
               J = J + 1
               GO TO 2450

            ELSE
               QFLAG = .FALSE.
            END IF

         ELSE IF (CH.EQ.CHAR(0)) THEN
            LXGTBP = .FALSE.
            RETURN

         END IF

      END IF

 2450 GO TO 2440

 2460 CONTINUE
      JLINE = J + 1
      LXGTBP = .TRUE.
      RETURN

      END
