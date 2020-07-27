C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INFACT (ML, MS, IFOUND, FACT, IIN, N19, N20, FACTOR,
     &   NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
C***********************************************************************

C  SUBROUTINE INFACT = ENTERS FACTORS INTO THE DATABASE

C***********************************************************************

      DIMENSION FACTOR (ML), IIN (IFOUND)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION LINKL (2, ML), LINKS (2, MS)

      LOGICAL ADDLNK

      IF (FACT .LE. 0.)FACT = 1.0
      DO 110 I = 1, IFOUND
         JJ = IIN (I)

C  APPLY INTERVALS TO A SIDE

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

C  INPUT A FACTOR FOR A LINE

            CALL LTSORT (ML, LINKL, JJ, IPNTR, ADDLNK)
            IF ( (JJ .GT. N19) .OR. (IPNTR .LE. 0)) THEN
               WRITE ( * , 10010)JJ
            ELSE
               FACTOR (IPNTR) = FACT
            ENDIF
         ENDIF
  110 CONTINUE

      RETURN

10000 FORMAT (' SIDE NO:', I5, ' IS NOT IN THE DATABASE',  / ,
     &   ' THUS NO INTERVALS CAN BE ENTERED')
10010 FORMAT (' LINE NO:', I5, ' IS NOT IN THE DATABASE',  / ,
     &   ' THUS NO INTERVAL CAN BE ENTERED')

      END
