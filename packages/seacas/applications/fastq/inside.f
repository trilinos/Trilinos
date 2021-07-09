C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INSIDE (MS, N3, N4, N20, JJ, IIN, IFOUND, ISIDE, NLPS,
     &   IFLINE, ILLIST, LINKS, NHOLDS, IHOLDS, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE INSIDE = INPUTS A SIDE INTO THE DATABASE

C***********************************************************************

      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (3 * MS)
      DIMENSION LINKS (2, MS)
      DIMENSION IIN (IFOUND), IHOLDS (2, MS)

      LOGICAL MERGE, NOROOM, ADDLNK

      IZ = 0
      NOROOM = .TRUE.

C  ZERO OUT THE LINK ARRAY IF NEEDED

      IF (JJ .GT. N20) THEN
         N20 = JJ

C  FIND THE CORRECT LINE NUMBER IF MERGING

      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (MS, LINKS, JJ, IPNTR, ADDLNK)
         JHOLD = JJ
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDS)NHOLDS = JHOLD
            CALL LTSORT (MS, IHOLDS, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N20 + 1
               N20 = N20 + 1
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MS, IHOLDS, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF

C  INPUT THE SIDE DATA INTO THE DATABASE

      N3 = N3 + 1
      J = N3
      IF (J .GT. MS)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MS, LINKS, JJ, J, ADDLNK)
      ISIDE (J) = JJ
      IFLINE (J) = N4 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 110
         N4 = N4 + 1
         IF (N4 .GT. MS * 3)RETURN
         ILLIST (N4) = JJ
  100 CONTINUE
  110 CONTINUE
      NLPS (J) = N4 - IFLINE (J) + 1
      IF (NLPS (J) .LT. 1) THEN
         WRITE ( * , 10010)J
         CALL LTSORT (MS, LINKS, ISIDE (J), IZ, ADDLNK)
      ENDIF

      NOROOM = .FALSE.
      RETURN

10000 FORMAT ('   OLD SIDE NO:', I5, '   TO NEW SIDE NO:', I5)
10010 FORMAT (' SIDE:', I5, ' HAS LESS THAN ONE LINE', / ,
     &   ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
      END
