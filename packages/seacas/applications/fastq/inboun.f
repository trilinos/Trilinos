C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INBOUN (MDIM, JJ, IFOUND, IIN, N1, N2, N3, N2OLD,
     &   MERGE, NOROOM, NEWNUM, NHOLD, IHOLD, IFLAG, INUM, IFIRST,
     &   LIST, LINK, IWT, JHOLD, ADDOLD)
C***********************************************************************

C  SUBROUTINE INBOUN = INPUTS AND LINKS BOUNDARY FLAG INFORMATION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     READ = READS AND/OR MERGES QMESH CARD FILE(S)
C   LINKBC = LINKS BOUNDARY FLAGS TO ENTITIES

C***********************************************************************

      DIMENSION IHOLD(2, MDIM), IFLAG(MDIM)
      DIMENSION INUM(MDIM), IFIRST(MDIM)
      DIMENSION LIST(2, MDIM), LINK(2, MDIM)
      DIMENSION IIN(IFOUND), IWT(3, MDIM)

      LOGICAL NOROOM, MERGE, NEWNUM, ADDOLD, ADDLNK

      NEWNUM = .FALSE.
      NOROOM = .FALSE.
      IZ = 0
      IOLD = 0

      IF ((JJ.LE.0) .OR. (JJ .GT. 10000)) THEN
         WRITE(*, 10000) JJ
         RETURN
      END IF

      NOROOM = .TRUE.

C  GET THE CORRECT FLAG ID

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

C  GET THE OLD LOCATION OF THE FLAG IF IT IS THERE

      ADDLNK = .FALSE.
      CALL LTSORT (MDIM, LINK, JJ, IOLD, ADDLNK)

C  IF THE FLAG CURRENTLY EXISTS,  SHIFT THE FLAG DATA TO THE END
C  OF THE CURRENT FLAG LIST,  TO FACILITATE ADDING MORE ENTRIES
C  TO THE OLD FLAG

      IF (IOLD .GT. 0) THEN

C  SHIFT THE OLD DEFINITION TO THE END OF THE LIST

         IFLAG(N2 + 1) = IFLAG(IOLD)
         INUM(N2 + 1) = INUM(IOLD)
         IWT(1, N2 + 1) = IWT(1, IOLD)
         IWT(2, N2 + 1) = IWT(2, IOLD)
         IWT(3, N2 + 1) = IWT(3, IOLD)
         IFIRST(N2 + 1) = N3 + 1

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

C  SLIDE ALL TRAILING FLAGS OVER TO FILL THE GAP IN THE LIST
C  RESORT THE POINTER ARRAY AS THE LIST FILLS AND STEP N2OLD
C  DOWN A NOTCH SO THESE WILL BE RELINKED IF NECESSARY

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

         N3 = IFIRST(N2) + INUM(N2) - 1
         DO 120 I = NLIST1, N3
            IKOUNT = I + KOUNT
            LIST(1, I) = LIST(1, IKOUNT)
            LIST(2, I) = LIST(2, IKOUNT)
  120    CONTINUE
      ELSE

C  LINK UP THE FLAG IN THE NEXT OPEN SLOT

         ADDLNK = .TRUE.
         N2 = N2 + 1
         CALL LTSORT (MDIM, LINK, JJ, N2, ADDLNK)
         IFLAG(N2) = JJ
         IFIRST(N2) = N3 + 1
         IWT(1, N2) = 0
         IWT(2, N2) = 0
         IWT(3, N2) = 0
      END IF

C  READ IN THE NEW FLAG LIST

      DO 150 I = 1, IFOUND
         JJ = IIN(I)
         IF (JJ .EQ. 0) GO TO 160

C  CHECK TO MAKE SURE THIS ENTITY IS NOT ALREADY IN THE LIST

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

      INUM(N2) = N3 - IFIRST(N2) + 1
      NOROOM = .FALSE.
      RETURN

10000 FORMAT (' A FLAG NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS FLAG WILL NOT BE INPUT INTO DATABASE')
      END
