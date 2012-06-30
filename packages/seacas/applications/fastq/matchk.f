C $Id: matchk.f,v 1.1 1990/11/30 11:11:57 gdsjaar Exp $
C $Log: matchk.f,v $
C Revision 1.1  1990/11/30 11:11:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]MATCHK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION MATCHK (MXND, I1, I2, J1, J2, LXN)
C***********************************************************************
C
C  FUNCTION MATCHK = CHECKS THE CURRENT COLAPSED LINES TO SEE IF THEY
C                    CAN BE JOINED WITHOUT AFFECTING THE BOUNDARY.
C                    I1 & I2 MAY END UP SWITCHED WITH J1 & J2.
C
C***********************************************************************
C
      DIMENSION LXN (4, MXND)
C
      IF ( (LXN (2, I1) .LT. 0) .OR. (LXN (2, I2) .LT. 0) .OR.
     &   (LXN (2, J1) .LT. 0) .OR. (LXN (2, J2) .LT. 0) ) THEN
C
C  FIRST CHECK FOR COMPLETELY HOOKED BOUNDARY LINES.
C
         IF ((LXN (2, J1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) THEN
            MATCHK = .FALSE.
         ELSEIF ( ((LXN (2, I1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) .OR.
     &      ((LXN (2, I2) .LT. 0) .AND. (LXN (2, J1) .LT. 0)))
     &      THEN
            MATCHK = .FALSE.
         ELSE
            MATCHK = .TRUE.
         ENDIF
      ELSE
         MATCHK = .TRUE.
      ENDIF
C
      RETURN
C
      END
