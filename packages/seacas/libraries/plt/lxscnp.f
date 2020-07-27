C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXSCNP(DELIM,STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) STR
      INTEGER NS
      CHARACTER CH
      LOGICAL QFLAG
      CHARACTER QCH

      NS = 0
      J = JLINE
      LP = 0
      QFLAG = .FALSE.
 2600 CONTINUE
      CH = ILINE(J:J)
      IF (LP.EQ.0) THEN
         ID = INDEX(DELIM,CH)
         IF (ID.GT.0) THEN
            GO TO 2620

         END IF

      END IF

      IF (.NOT.QFLAG) THEN
         IF (LP.EQ.0) THEN
            ID = INDEX(DELIM,CH)
            IF (ID.GT.0) THEN
               GO TO 2620

            END IF

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.'(' .OR. CH.EQ.'[' .OR. CH.EQ.'{') THEN
            LP = LP + 1

         ELSE IF (CH.EQ.')' .OR. CH.EQ.']' .OR. CH.EQ.'}') THEN
            LP = LP - 1

         ELSE IF (CH.EQ.'''' .OR. CH.EQ.'"' .OR. CH.EQ.CHAR(96)) THEN
            QFLAG = .TRUE.
            QCH = CH
         END IF

         IF (LP.LT.0) THEN
            GO TO 2620

         END IF

         NS = NS + 1
         STR(NS:NS) = CH

      ELSE
         NS = NS + 1
         STR(NS:NS) = CH
         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.QCH) THEN
            IF (CH.EQ.ILINE(J+1:J+1)) THEN
               NS = NS + 1
               STR(NS:NS) = CH
               J = J + 2
               GO TO 2610

            ELSE
               QFLAG = .FALSE.
            END IF

         END IF

      END IF

      J = J + 1
 2610 GO TO 2600

 2620 CONTINUE
      JLINE = J
      LXSCNP = ((LP.EQ.0) .AND. .NOT.QFLAG)
      RETURN

      END
