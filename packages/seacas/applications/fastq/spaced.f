C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES,
     &   ICOMB, ITEST, LTEST, ERR)
C***********************************************************************

C  SUBROUTINE SPACED = COUNTS THE INTERVAL SPACINGS FOR A COMBINATION

C***********************************************************************

      DIMENSION ICOMB (MXCORN), LCORN (MXCORN)
      DIMENSION ITEST (ILEN), LTEST (ILEN)
      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR

      ERR = .TRUE.
      KLEN = 0
      KOUNTC = 0

      DO 100 I = 1, NCORN

         IF (ICOMB (I) .EQ. 1) THEN
            KLEN = KLEN + 1
            IF (KLEN .GT. ILEN) THEN
               CALL MESAGE ('PROBLEMS IN SPACED - COUNTERS DON''T '//
     &            'MATCH DATA')
               RETURN
            ENDIF

            ITEST (KLEN) = LCORN(I)
            LTEST (KLEN) = KOUNTC
            KOUNTC = LNODES (7, LCORN(I))
         ELSE
            KOUNTC = KOUNTC + LNODES (7, LCORN (I))
         ENDIF
  100 CONTINUE

C  NOW ADD THE REMAINING KOUNTC ONTO THE FRONT

      LTEST (1) = LTEST (1) + KOUNTC

C  NOW SWITCH THE COUNTS TO BE FOLLOWING THE CORNERS INSTEAD OF
C  BEFORE THE CORNERS

      IHOLD = LTEST (1)
      DO 110 I = 2, KLEN
         LTEST (I - 1) = LTEST (I)
  110 CONTINUE
      LTEST (KLEN) = IHOLD

      ERR = .FALSE.
      RETURN

      END
