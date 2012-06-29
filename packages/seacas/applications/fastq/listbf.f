C $Id: listbf.f,v 1.4 1998/07/14 18:19:20 gdsjaar Exp $
C $Log: listbf.f,v $
C Revision 1.4  1998/07/14 18:19:20  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.3  1998/07/14 17:42:17  gdsjaar
C *** empty log message ***
C
C Revision 1.2  1998/04/16 05:06:44  gdsjaar
C Changed "X" to "1X" in format statement
C
C Revision 1.1.1.1  1990/11/30 11:11:20  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:11:19  gdsjaar
c Initial revision
c 
CC* FILE: [.MAIN]LISTBF.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LISTBF (MDIM, N, CHOICE, LINK, IFLAG, INUM, IFIRST,
     &   LIST, IWT)
C***********************************************************************
C
C  SUBROUTINE LISTBF = LISTS BOUNDARY CONDITIONS BY FLAG NUMBERS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     LIST = LISTS POINTS,  LINES,  REGIONS,  SCHEMES,  AND BOUNDARY
C            DEFINITIONS
C
C***********************************************************************
C
      DIMENSION LINK (2, MDIM), IFLAG (MDIM), INUM (MDIM), IFIRST (MDIM)
      DIMENSION LIST (2, MDIM), IWT (3, MDIM)
C
      CHARACTER CHOICE*7
C
      LOGICAL ADDLNK, EXTRA, FOUND
C
      ADDLNK = .FALSE.
      FOUND = .FALSE.
C
      IF (CHOICE (1:5) .EQ. 'POINT') THEN
         WRITE (*, 10000)
      ENDIF
      IF (N .GT. 0) THEN
         DO 110 I = 1, N
            CALL LTSORT (MDIM, LINK, I, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               FOUND = .TRUE.
               EXTRA = .FALSE.
               L1 = IFIRST (IPNTR)
               L2 = L1 + INUM (IPNTR) - 1
  100          CONTINUE
               IF ( (L2 - L1 + 1) .GT. 6) THEN
                  L2 = L1 + 5
                  IF (EXTRA) THEN
                     WRITE (*, 10010) (LIST (1, L), L = L1, L2)
                  ELSE
                     IF (IWT (1, IPNTR) .EQ. 0) THEN
                        WRITE (*, 10020) IFLAG (IPNTR), CHOICE,
     &                     (LIST (1, L),  L = L1, L2)
                     ELSE
                        IF (IWT (3, IPNTR) .EQ. 0) THEN
                           WRITE (*, 10030) IFLAG (IPNTR), CHOICE,
     &                        IWT (1, IPNTR),  IWT (2, IPNTR),
     &                        (LIST (1, L), L = L1, L2)
                        ELSE
                           WRITE (*, 10040) IFLAG (IPNTR), CHOICE,
     &                        IWT (1, IPNTR),  IWT (2, IPNTR),
     &                        IWT (3, IPNTR),  (LIST (1, L), L = L1, L2)
                        ENDIF
                     ENDIF
                     EXTRA = .TRUE.
                  ENDIF
                  L1 = L2 + 1
                  L2 = IFIRST (IPNTR) + INUM (IPNTR)-1
                  GOTO 100
               ELSE
                  IF (EXTRA) THEN
                     WRITE (*, 10010) (LIST (1, L), L = L1, L2)
                  ELSE
                     IF (IWT (1, IPNTR) .EQ. 0) THEN
                        WRITE (*, 10020) IFLAG (IPNTR), CHOICE,
     &                     (LIST (1, L), L = L1, L2)
                     ELSE
                        IF (IWT (3, IPNTR) .EQ. 0) THEN
                           WRITE (*, 10030) IFLAG (IPNTR), CHOICE,
     &                        IWT (1, IPNTR), IWT (2, IPNTR),
     &                        (LIST (1, L), L = L1, L2)
                        ELSE
                           WRITE (*, 10040) IFLAG (IPNTR), CHOICE,
     &                        IWT (1, IPNTR), IWT (2, IPNTR),
     &                        IWT (3, IPNTR),  (LIST (1, L), L = L1, L2)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  110    CONTINUE
      ELSE
         FOUND = .TRUE.
         WRITE (*, 10050) CHOICE
      ENDIF
      IF ( .NOT. FOUND) THEN
         WRITE (*, 10050) CHOICE
      ENDIF
      RETURN
C
10000 FORMAT ('  FLAG    BOUN.  WEIGHT FIRST  FIRST', /,
     &   ' NUMBER   TYPE    SIDE  WT PNT WT LIN     POINT OR LINE '//
     &   'LISTING', /, ' ------ -------- ------ ------ ------ ',
     &   '----------------------------------')
10010 FORMAT (' ', 37X, 6I6)
10020 FORMAT (' ',1X, I5, 2X, A7, 2X, '-----', 2X, '-----', 2X,
     &   '-----', 2X, 6I6)
10030 FORMAT (' ',1X, I5, 2X, A7, 2X, I5, 2X, I5, 2X, '-----', 2X, 6I6)
10040 FORMAT (' ',1X, I5, 2X, A7, 2X, I5, 2X, I5, 2X, I5, 2X, 6I6)
10050 FORMAT (' *** NO ', A7, ' FLAGS IN THE CURRENT DATABASE ***')
      END
