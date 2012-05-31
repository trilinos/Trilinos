C $Id: fixlxn.f,v 1.1 1990/11/30 11:07:27 gdsjaar Exp $
C $Log: fixlxn.f,v $
C Revision 1.1  1990/11/30 11:07:27  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]FIXLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNNOLD, LLLOLD, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE FIXLXN = FIXES THE ADDITIONS TO LXN
C
C***********************************************************************
C
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, NOROOM
C
C     RE-SETUP AVAILABLE LXN-SPACE LINKS
C
      IOLD = 0
      NAVAIL = 0
      DO 100 I = 1, NNNOLD
         IF (LXN (1, I).EQ.0)THEN
            IF (IOLD.LE.0)THEN
               IAVAIL = I
            ELSE
               LXN (4, IOLD) = I
            ENDIF
            IOLD = I
            NAVAIL = NAVAIL + 1
         ENDIF
  100 CONTINUE
      IF (IOLD.LE.0)THEN
         IAVAIL = NNN + 1
      ELSE
         LXN (4, IOLD) = NNN + 1
      ENDIF
      NAVAIL = NAVAIL + (MXND - NNN)
      IF (NNN.LT.MXND)THEN
         NNN1 = NNN + 1
         DO 110 I = NNN1, MXND
            LXN (1, I) = 0
            LXN (2, I) = 0
            LXN (3, I) = 0
            LXN (4, I) = I + 1
  110    CONTINUE
      ENDIF
C
C     COMPLETE LXN ARRAYS FOR ANY NEW LINES
C
      DO 130 L = LLLOLD + 1, LLL
         DO 120 I = 1, 2
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NXL (I, L),
     &         L, NNN, ERR, NOROOM)
            IF (ERR)THEN
               CALL MESAGE ('ERROR IN FIXLXN - NXL TABLE GENERATION')
               GOTO 140
            ELSEIF (NOROOM) THEN
               GOTO 140
            ENDIF
  120    CONTINUE
  130 CONTINUE
C
  140 CONTINUE
      RETURN
C
      END
