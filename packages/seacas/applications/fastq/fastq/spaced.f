C $Id: spaced.f,v 1.1 1990/11/30 11:16:18 gdsjaar Exp $
C $Log: spaced.f,v $
C Revision 1.1  1990/11/30 11:16:18  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SPACED.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES,
     &   ICOMB, ITEST, LTEST, ERR)
C***********************************************************************
C
C  SUBROUTINE SPACED = COUNTS THE INTERVAL SPACINGS FOR A COMBINATION
C
C***********************************************************************
C
      DIMENSION ICOMB (MXCORN), LCORN (MXCORN)
      DIMENSION ITEST (ILEN), LTEST (ILEN)
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR
C
      ERR = .TRUE.
      KLEN = 0
      KOUNTC = 0
C
      DO 100 I = 1, NCORN
C
         IF (ICOMB (I) .EQ. 1) THEN
            KLEN = KLEN + 1
            IF (KLEN .GT. ILEN) THEN
               CALL MESAGE ('PROBLEMS IN SPACED - COUNTERS DON''T '//
     &            'MATCH DATA')
               RETURN
            ENDIF
C
            ITEST (KLEN) = LCORN(I)
            LTEST (KLEN) = KOUNTC
            KOUNTC = LNODES (7, LCORN(I))
         ELSE
            KOUNTC = KOUNTC + LNODES (7, LCORN (I))
         ENDIF
  100 CONTINUE
C
C  NOW ADD THE REMAINING KOUNTC ONTO THE FRONT
C
      LTEST (1) = LTEST (1) + KOUNTC
C
C  NOW SWITCH THE COUNTS TO BE FOLLOWING THE CORNERS INSTEAD OF
C  BEFORE THE CORNERS
C
      IHOLD = LTEST (1)
      DO 110 I = 2, KLEN
         LTEST (I - 1) = LTEST (I)
  110 CONTINUE
      LTEST (KLEN) = IHOLD
C
      ERR = .FALSE.
      RETURN
C
      END
