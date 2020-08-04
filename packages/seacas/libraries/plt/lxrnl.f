C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXRNL(VAL,N,CH)
      DOUBLE PRECISION VAL(*)
      CHARACTER CH
      LOGICAL LDUM,LXGTCH,LXGTWH,LXREAL
      DOUBLE PRECISION XX,XY

      I = 1
 2320 IF (.NOT. (I.LE.N)) GO TO 2340
      VAL(I) = 0.
      I = I + 1
      GO TO 2320

 2340 CONTINUE
      N = 0
      LXRNL = .TRUE.
 2350 CONTINUE
      LDUM = LXGTWH(CH)
      IF (LXREAL(VAL(N+1),CH)) THEN
         N = N + 1
         LDUM = LXGTWH(CH)
         IF (LXGTCH('#',CH)) THEN
            RETURN

         END IF

         IF (LXGTCH(',',CH)) THEN
            GO TO 2360

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            RETURN

         END IF

         IF (LXGTCH('*',CH)) THEN
            LDUM = LXGTWH(CH)
            XX = VAL(N)
            IF (.NOT.LXREAL(XY,CH)) THEN
               LXRNL = .FALSE.
               RETURN

            END IF

            M = INT(XX + .1)
            N0 = N
 2380       IF (.NOT. (N.LT.M+N0)) GO TO 2400
            VAL(N) = XY
            N = N + 1
            GO TO 2380

 2400       CONTINUE
            LDUM = LXGTWH(CH)
            N = N0 + MAX(M-1,0)
            IF (LXGTCH(',',CH)) THEN
               GO TO 2360

            END IF

            IF (LXGTCH('#',CH)) THEN
               RETURN

            END IF

            IF (CH.EQ.CHAR(0)) THEN
               RETURN

            END IF

         END IF

      ELSE IF (LXGTCH(',',CH)) THEN
         VAL(N+1) = 0.
         N = N + 1

      ELSE IF (LXGTCH('#',CH)) THEN
         RETURN

      ELSE IF (CH.EQ.CHAR(0)) THEN
         IF (N.EQ.0) THEN
            RETURN

         END IF

         N = N + 1
         VAL(N) = 0.
         RETURN

      ELSE
         LXRNL = .FALSE.
         RETURN

      END IF

 2360 GO TO 2350

      END
