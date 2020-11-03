C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ININTR (ML, MS, IFOUND, NEWINT, IIN, N19, N20, NINT,
     &   NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
C***********************************************************************

C  SUBROUTINE ININTR = ENTERS INTERVALS INTO THE DATABASE

C***********************************************************************

      DIMENSION NINT (ML), IIN (IFOUND)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION LINKL (2, ML), LINKS (2, MS)

      LOGICAL ADDLNK

      DO 110 I = 1, IFOUND
         JJ = IIN (I)

C  APPLY INTERVALS TO A SIDE

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

C  INPUT INTERVALS ON A LINE

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

      RETURN

10000 FORMAT (' SIDE NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO INTERVALS CAN BE ENTERED')
10010 FORMAT (' LINES IN SIDE NO:', I5,
     &   ' CANNOT HAVE INTERVALS OF:', I8)
10020 FORMAT (' LINE NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO INTERVAL CAN BE ENTERED')
10030 FORMAT (' LINE NO:', I5, ' CANNOT HAVE AN INTERVAL OF:', I8)

      END
