C $Id: inintr.f,v 1.1 1990/11/30 11:09:46 gdsjaar Exp $
C $Log: inintr.f,v $
C Revision 1.1  1990/11/30 11:09:46  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ININTR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ININTR (ML, MS, IFOUND, NEWINT, IIN, N19, N20, NINT,
     &   NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
C***********************************************************************
C
C  SUBROUTINE ININTR = ENTERS INTERVALS INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION NINT (ML), IIN (IFOUND)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION LINKL (2, ML), LINKS (2, MS)
C
      LOGICAL ADDLNK
C
      DO 110 I = 1, IFOUND
         JJ = IIN (I)
C
C  APPLY INTERVALS TO A SIDE
C
         IF (JJ.LT.0) THEN
            JJ = - JJ
            CALL LTSORT (MS, LINKS, JJ, IPNTR, ADDLNK)
            IF ( (JJ.GT.N20) .OR. (IPNTR.LE.0)) THEN
               WRITE (*, 10000) JJ
            ELSEIF (NEWINT.LT.0) THEN
               WRITE (*, 10010) JJ, IIN (2)
            ELSE
               DO 100 JJJ = IFLINE (IPNTR), IFLINE (IPNTR) +
     &            NLPS (IPNTR)-1
                  CALL LTSORT (ML, LINKL, ILLIST (JJJ), JPNTR,
     &               ADDLNK)
                  IF ( (ILLIST (JJJ) .GT. N19) .OR.
     &               (IPNTR.LE.0)) THEN
                     WRITE (*, 10020) ILLIST (JJJ)
                  ELSE
                     NINT (JPNTR) = NEWINT
                  ENDIF
  100          CONTINUE
            ENDIF
         ELSE
C
C  INPUT INTERVALS ON A LINE
C
            CALL LTSORT (ML, LINKL, JJ, IPNTR, ADDLNK)
            IF ( (JJ.GT.N19) .OR. (IPNTR.LE.0)) THEN
               WRITE (*, 10020) JJ
            ELSEIF (NEWINT .LT. 0) THEN
               WRITE (*, 10030) JJ, NEWINT
            ELSE
               NINT (IPNTR) = NEWINT
            ENDIF
         ENDIF
  110 CONTINUE
C
      RETURN
C
10000 FORMAT (' SIDE NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO INTERVALS CAN BE ENTERED')
10010 FORMAT (' LINES IN SIDE NO:', I5,
     &   ' CANNOT HAVE INTERVALS OF:', I8)
10020 FORMAT (' LINE NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO INTERVAL CAN BE ENTERED')
10030 FORMAT (' LINE NO:', I5, ' CANNOT HAVE AN INTERVAL OF:', I8)
C
      END
