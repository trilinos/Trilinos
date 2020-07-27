C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LINKBC (MDIM, MS, I1, I2, NBOUN, N1, N2, N3, N20,
     &   IFLAG, IFLIST, NEPS, LIST, NLPS, IFLINE, ILLIST, IBOUN,
     &   LINKF, IWT, LINKE, LINKS, SIDEOK, NOROOM)
C***********************************************************************

C  SUBROUTINE LINKBC = LINKS UP ALL BOUNDARY FLAG LISTS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     READ = READS AND/OR MERGES FASTQ FILE(S)

C***********************************************************************

C  VARIABLES USED:
C     I1     = THE FIRST FLAG TO BE LINKED
C     I2     = THE LAST FLAG TO BE LINKED
C     IFLAG  = THE ARRAY OF FLAGS
C     IFLIST = THE FIRST ENTITY IN LIST TO BE ASSOCIATED WITH A FLAG
C     NEPS   = THE NUMBER OF ENTITIES IN LIST THAT GO WITH A FLAG
C     LIST   = THE LIST OF ENTITIES
C     LINK   = THE LINK TO THE FLAG LIST
C     IBOUN  = THE LINK FROM THE ENTITY TO THE FLAGS
C     MDIM   = THE DIMENSIONING PARAMETER FOR THE LIST
C     SIDEOK = .FALSE. IF IT IS NOT POSSIBLE TO EXPAND SIDES  (POINBC'S)

C***********************************************************************

      DIMENSION IFLAG(MDIM), IFLIST(MDIM), NEPS(MDIM), LIST(2, MDIM)
      DIMENSION LINKF(2, MDIM), IBOUN(MDIM)
      DIMENSION LINKE(2, MDIM), LINKS(2, MS)
      DIMENSION NLPS(MS), IFLINE(MS), ILLIST(MS*3), IWT(3, MDIM)

      LOGICAL ADDLNK, MERGE, ADDOLD, NOROOM, SIDEOK, NEWNUM

      IZ = 0
      ADDLNK = .FALSE.
      MERGE = .FALSE.
      NOROOM = .FALSE.

      IF (SIDEOK) THEN

C  EXPAND ALL THE SIDES (SETS) TO THEIR RESPECT LINES (ENTITIES)

         DO 130 I = I1, I2
  100       CONTINUE
            CALL LTSORT (MDIM, LINKF, IFLAG(I), II, ADDLNK)
            IF (II .GT. 0) THEN

C  THE FLAG HAS BEEN FOUND

               IFLAG1 = IFLAG(II)
               J1 = IFLIST(II)
               J2 = J1 + NEPS(II) - 1
               DO 120 J = J1, J2
                  JJ = LIST(1, J)
                  IF (JJ .LT. 0) THEN

C REMOVE THE SIDE FROM THE FLAG LIST

                     NEPS(II) = NEPS(II) - 1
                     DO 110 K = J, J2 - 1
                        LIST(1, K) = LIST(1, K + 1)
                        LIST(2, K) = LIST(2, K + 1)
  110                CONTINUE

C  IF THE SIDE EXISTS,  REPLACE IT WITH THE LINES IT REPRESENTS

                     JJ = -JJ
                     CALL LTSORT (MS, LINKS, JJ, IPNTR, ADDLNK)
                     IF ((JJ .GT. N20) .OR. (IPNTR .LE. 0)) THEN
                        WRITE(*, 10000) JJ
                     ELSE
                        CALL INBOUN (MDIM, IFLAG1, NLPS(IPNTR),
     &                     ILLIST(IFLINE(IPNTR)), N1, N2, N3, N2OLD,
     &                     MERGE, NOROOM, NEWNUM, IZ, LINKF, IFLAG,
     &                     NEPS, IFLIST, LIST, LINKF, IWT, IZ, ADDOLD)
                        IF (NOROOM) RETURN

C  NOW,  SEE IF THERE ARE ANY SIDES IN THE NEW I'TH FLAG'S LIST
C  NOTE THAT THE ONE FIXED HAS NOW BEEN ROTATED TO THE END OF THE LIST.

                        GOTO 100
                     ENDIF
                  ENDIF
  120          CONTINUE
            ENDIF
  130    CONTINUE
      ENDIF

C  ALL POSSIBLE SIDE EXPANSION HAS OCCURRED
C  NOW LINK UP ALL THE LINES

      DO 160 I = I1, I2
         CALL LTSORT (MDIM, LINKF, IFLAG(I), II, ADDLNK)
         IF (II .GT. 0) THEN

C  THE FLAG HAS BEEN FOUND

            IFLAG1 = IFLAG(II)
            J1 = IFLIST(II)
            J2 = J1 + NEPS(II) - 1
            DO 140 J = J1, J2
               JJ = LIST(1, J)
               IF (JJ .GT. 0) THEN
                  CALL LINKEN (MDIM, JJ, IFLAG1, IFLAG, IFLIST, NEPS,
     &               LIST, LINKF, LINKE, IBOUN, ADDLNK)
               ELSE
                  CALL MESAGE ('PROBLEMS ELIMINATING SIDES IN LINKBC')
               ENDIF
  140       CONTINUE
         ELSE
            DO 150 J = 1, NBOUN
               IF (IBOUN(J) .EQ. IFLAG(I)) IBOUN(J) = 0
  150       CONTINUE
         ENDIF
  160 CONTINUE
      RETURN

10000 FORMAT (' SIDE NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO BOUNDARY FLAGS CAN BE ENTERED ALONG THIS SIDE')
      END
