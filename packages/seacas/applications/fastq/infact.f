C $Id: infact.f,v 1.1 1990/11/30 11:09:35 gdsjaar Exp $
C $Log: infact.f,v $
C Revision 1.1  1990/11/30 11:09:35  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INFACT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INFACT (ML, MS, IFOUND, FACT, IIN, N19, N20, FACTOR,
     &   NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
C***********************************************************************
C
C  SUBROUTINE INFACT = ENTERS FACTORS INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION FACTOR (ML), IIN (IFOUND)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION LINKL (2, ML), LINKS (2, MS)
C
      LOGICAL ADDLNK
C
      IF (FACT .LE. 0.)FACT = 1.0
      DO 110 I = 1, IFOUND
         JJ = IIN (I)
C
C  APPLY INTERVALS TO A SIDE
C
         IF (JJ .LT. 0) THEN
            JJ = - JJ
            CALL LTSORT (MS, LINKS, JJ, IPNTR, ADDLNK)
            IF ( (JJ .GT. N20) .OR. (IPNTR .LE. 0)) THEN
               WRITE ( * , 10000)JJ
            ELSE
               DO 100 JJJ = IFLINE (IPNTR),
     &            IFLINE (IPNTR) + NLPS (IPNTR) - 1
                  CALL LTSORT (ML, LINKL, ILLIST (JJJ), JPNTR, ADDLNK)
                  IF ( (ILLIST (JJJ) .GT. N19) .OR. (IPNTR .LE. 0))
     &               THEN
                     WRITE ( * , 10010)ILLIST (JJJ)
                  ELSE
                     FACTOR (JPNTR) = FACT
                  ENDIF
  100          CONTINUE
            ENDIF
         ELSE
C
C  INPUT A FACTOR FOR A LINE
C
            CALL LTSORT (ML, LINKL, JJ, IPNTR, ADDLNK)
            IF ( (JJ .GT. N19) .OR. (IPNTR .LE. 0)) THEN
               WRITE ( * , 10010)JJ
            ELSE
               FACTOR (IPNTR) = FACT
            ENDIF
         ENDIF
  110 CONTINUE
C
      RETURN
C
10000 FORMAT (' SIDE NO:', I5, ' IS NOT IN THE DATABASE',  / ,
     &   ' THUS NO INTERVALS CAN BE ENTERED')
10010 FORMAT (' LINE NO:', I5, ' IS NOT IN THE DATABASE',  / ,
     &   ' THUS NO INTERVAL CAN BE ENTERED')
C
      END
