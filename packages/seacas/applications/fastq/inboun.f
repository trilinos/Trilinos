C $Id: inboun.f,v 1.3 2004/01/21 05:18:40 gdsjaar Exp $
C $Log: inboun.f,v $
C Revision 1.3  2004/01/21 05:18:40  gdsjaar
C Initialized several variables identified by valgrind.
C
C Revision 1.2  1998/07/14 18:19:11  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:09:21  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:09:19  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]INBOUN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INBOUN (MDIM, JJ, IFOUND, IIN, N1, N2, N3, N2OLD,
     &   MERGE, NOROOM, NEWNUM, NHOLD, IHOLD, IFLAG, INUM, IFIRST,
     &   LIST, LINK, IWT, JHOLD, ADDOLD)
C***********************************************************************
C
C  SUBROUTINE INBOUN = INPUTS AND LINKS BOUNDARY FLAG INFORMATION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     READ = READS AND/OR MERGES QMESH CARD FILE(S)
C   LINKBC = LINKS BOUNDARY FLAGS TO ENTITIES
C
C***********************************************************************
C
      DIMENSION IHOLD(2, MDIM), IFLAG(MDIM)
      DIMENSION INUM(MDIM), IFIRST(MDIM)
      DIMENSION LIST(2, MDIM), LINK(2, MDIM)
      DIMENSION IIN(IFOUND), IWT(3, MDIM)
C
      LOGICAL NOROOM, MERGE, NEWNUM, ADDOLD, ADDLNK
C
      NEWNUM = .FALSE.
      NOROOM = .FALSE.
      IZ = 0
      IOLD = 0
C
      IF ((JJ.LE.0) .OR. (JJ .GT. 10000)) THEN
         WRITE(*, 10000) JJ
         RETURN
      END IF
C
      NOROOM = .TRUE.
C
C  GET THE CORRECT FLAG ID
C
      IF (JJ .GT. N1) THEN
         N1 = JJ
      ELSE IF (MERGE) THEN
         JHOLD = JJ
         ADDLNK = .FALSE.
         CALL LTSORT (MDIM, LINK, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLD)NHOLD = JHOLD
            CALL LTSORT (MDIM, IHOLD, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N1 + 1
               N1 = JJ
               NEWNUM = .TRUE.
               ADDLNK = .TRUE.
               CALL LTSORT (MDIM, LINK, JJ, IZ, ADDLNK)
               CALL LTSORT (MDIM, IHOLD, JHOLD, JJ, ADDLNK)
            END IF
         END IF
      END IF
      IF (N2 + 1 .GT. MDIM) RETURN
C
C  GET THE OLD LOCATION OF THE FLAG IF IT IS THERE
C
      ADDLNK = .FALSE.
      CALL LTSORT (MDIM, LINK, JJ, IOLD, ADDLNK)
C
C  IF THE FLAG CURRENTLY EXISTS,  SHIFT THE FLAG DATA TO THE END
C  OF THE CURRENT FLAG LIST,  TO FACILITATE ADDING MORE ENTRIES
C  TO THE OLD FLAG
C
      IF (IOLD .GT. 0) THEN
C
C  SHIFT THE OLD DEFINITION TO THE END OF THE LIST
C
         IFLAG(N2 + 1) = IFLAG(IOLD)
         INUM(N2 + 1) = INUM(IOLD)
         IWT(1, N2 + 1) = IWT(1, IOLD)
         IWT(2, N2 + 1) = IWT(2, IOLD)
         IWT(3, N2 + 1) = IWT(3, IOLD)
         IFIRST(N2 + 1) = N3 + 1
C
         IF (IOLD .LT. N2) THEN
            KOUNT = IFIRST(IOLD + 1) - IFIRST(IOLD)
         ELSE IF (IOLD .EQ. N2) THEN
            KOUNT = INUM(IOLD)
         ELSE
            CALL MESAGE ('IN INBOUN, ERROR REORDERING FLAG LIST')
            RETURN
         END IF
         NLIST1 = IFIRST(IOLD)
         NPLACE = N3
         DO 100 I = NLIST1, NLIST1 + INUM(IOLD) - 1
            NPLACE = NPLACE + 1
            IF (NPLACE .GT. MDIM) RETURN
            LIST(1, NPLACE) = LIST(1, I)
            LIST(2, NPLACE) = LIST(2, I)
  100    CONTINUE
C
C  SLIDE ALL TRAILING FLAGS OVER TO FILL THE GAP IN THE LIST
C  RESORT THE POINTER ARRAY AS THE LIST FILLS AND STEP N2OLD
C  DOWN A NOTCH SO THESE WILL BE RELINKED IF NECESSARY
C
         IF (N2OLD .GT. 0) N2OLD = N2OLD - 1
         ADDLNK = .TRUE.
         DO 110 I = IOLD, N2
            IFLAG(I) = IFLAG(I + 1)
            INUM(I) = INUM(I + 1)
            IFIRST(I) = IFIRST(I + 1) - KOUNT
            IWT(1, I) = IWT(1, I + 1)
            IWT(2, I) = IWT(2, I + 1)
            IWT(3, I) = IWT(3, I + 1)
            CALL LTSORT (MDIM, LINK, IFLAG(I), I, ADDLNK)
  110    CONTINUE
C
         N3 = IFIRST(N2) + INUM(N2) - 1
         DO 120 I = NLIST1, N3
            IKOUNT = I + KOUNT
            LIST(1, I) = LIST(1, IKOUNT)
            LIST(2, I) = LIST(2, IKOUNT)
  120    CONTINUE
      ELSE
C
C  LINK UP THE FLAG IN THE NEXT OPEN SLOT
C
         ADDLNK = .TRUE.
         N2 = N2 + 1
         CALL LTSORT (MDIM, LINK, JJ, N2, ADDLNK)
         IFLAG(N2) = JJ
         IFIRST(N2) = N3 + 1
         IWT(1, N2) = 0
         IWT(2, N2) = 0
         IWT(3, N2) = 0
      END IF
C
C  READ IN THE NEW FLAG LIST
C
      DO 150 I = 1, IFOUND
         JJ = IIN(I)
         IF (JJ .EQ. 0) GO TO 160
C
C  CHECK TO MAKE SURE THIS ENTITY IS NOT ALREADY IN THE LIST
C
         DO 130 K = IFIRST(N2), N3
            IF (LIST(1, K) .EQ. JJ) GO TO 140
  130    CONTINUE
         N3 = N3 + 1
         IF (N3 .GT. MDIM) RETURN
         LIST(1, N3) = JJ
         LIST(2, N3) = 0
  140    CONTINUE
  150 CONTINUE
  160 CONTINUE
C
      INUM(N2) = N3 - IFIRST(N2) + 1
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT (' A FLAG NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS FLAG WILL NOT BE INPUT INTO DATABASE')
      END
