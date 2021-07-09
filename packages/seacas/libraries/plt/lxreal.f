C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXREAL(VALUE,CH)
      LOGICAL LXSET,LXNUMB
      LOGICAL PREF,POSF
      DOUBLE PRECISION VALUE,PREDIG,ESIGN
      DOUBLE PRECISION FN
      CHARACTER CH
      INTEGER ND

      LXREAL = .FALSE.
      ESIGN = 1.
      PREDIG = 0.
      ISIGN = 1
      ISAVE = LXSV()
      IF (LXSET('+-',CH)) THEN
         IF (CH.EQ.'-') THEN
            ISIGN = -1
         END IF

      END IF

      PREF = LXNUMB(PREDIG,ND,CH)
      POSF = .FALSE.
      IF (LXSET('.',CH)) THEN
         IF (LXNUMB(FN,ND,CH)) THEN
            POSF = .TRUE.
            PREDIG = PREDIG + FN*10.** (DBLE(-ND))
         END IF

      END IF

      PREDIG = PREDIG*ISIGN
      IF (.NOT. (PREF.OR.POSF)) THEN
         CALL LXRS(ISAVE)
         RETURN

      END IF

      IF (LXSET('EeDdQq',CH)) THEN
         IF (LXSET('+-',CH)) THEN
            IF (CH.EQ.'-') THEN
               ESIGN = -1.
            END IF

         END IF

         IF (LXNUMB(FN,ND,CH)) THEN
            PREDIG = PREDIG*10.** (ESIGN*FN)

         ELSE
            CALL LXRS(ISAVE)
            RETURN

         END IF

      END IF

      VALUE = PREDIG
      LXREAL = .TRUE.
      RETURN

      END
