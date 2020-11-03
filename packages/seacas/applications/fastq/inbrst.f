C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INBRST (MS, MR, N5, N6, N21, N23, JJ, IMTRL, JC, IIN,
     &   IFOUND, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST, LINKB,
     &   LINKM, NHOLDM, IHOLDM, NHOLDB, IHOLDB, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE INBRST = INPUTS A BAR SET INTO THE DATABASE

C***********************************************************************

      DIMENSION IBARST (MS), JMAT (MS), JCENT (MS), NLPB (MS)
      DIMENSION JFLINE (MS)
      DIMENSION JLLIST (3 * MS), LINKB (2, MS), LINKM (2, MS + MR)
      DIMENSION IHOLDM (2,  (MS + MR)), IHOLDB (2, MS)
      DIMENSION IIN (IFOUND)

      LOGICAL MERGE, NOROOM, ADDLNK

      IZ = 0
      NOROOM = .TRUE.
      N22 = 0

C  ZERO OUT THE LINK ARRAY IF NEEDED

      IF (JJ .GT. N21) THEN
         N21 = JJ

C  FIND THE CORRECT BAR SET NUMBER IF MERGING

C  SET UP POINTERS FOR MERGING DATA

      ELSEIF (MERGE) THEN
         JHOLD = JJ
         CALL LTSORT (MS, LINKB, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDB)NHOLDB = JHOLD
            CALL LTSORT (MS, IHOLDB, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N22 + 1
               N22 = JJ
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MS, IHOLDB, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF

C  INPUT THE BAR SET DATA INTO THE DATABASE

      N5 = N5 + 1
      J = N5
      IF (J .GT. MS)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MS, LINKB, JJ, J, ADDLNK)
      IBARST (J) = JJ
      JCENT (J) = JC
      JFLINE (J) = N6 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 110
         N6 = N6 + 1
         IF (N6 .GT. MS * 3)RETURN
         JLLIST (N6) = JJ
  100 CONTINUE
  110 CONTINUE
      NLPB (J) = N6 - JFLINE (J) + 1
      IF (NLPB (J) .LT. 1) THEN
         WRITE ( * , 10010)J
         CALL LTSORT (MS, LINKB, JJ, IZ, ADDLNK)
      ENDIF
      ADDLNK = .FALSE.

C  LINK UP THE MATERIAL

C  ZERO THE LINK ARRAY IF NEEDED

      IF (IMTRL .GT. N23) THEN
         N23 = IMTRL

C  SET UP POINTERS FOR MERGING DATA

      ELSEIF (MERGE) THEN
         JHOLD = IMTRL
         CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
         IF (IPNTR .NE. 0) THEN
            IF (JHOLD .GT. NHOLDM)NHOLDM = JHOLD
            ADDLNK = .FALSE.
            CALL LTSORT ( (MS + MR), IHOLDM, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               IMTRL = IPNTR
            ELSE
               IMTRL = N23 + 1
               N23 = IMTRL
               WRITE ( * , 10010)JHOLD, IMTRL
               ADDLNK = .TRUE.
               CALL LTSORT ( (MS + MR), IHOLDM, JHOLD, IMTRL, ADDLNK)
            ENDIF
         ENDIF
      ENDIF

C  ADD THE MATERIAL INTO THE DATABASE

      NOROOM = .FALSE.
      ADDLNK = .FALSE.
      CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
      ADDLNK = .TRUE.
      IF (IPNTR .GT. 0) THEN
         CALL MESAGE (' ')
         WRITE ( * , 10020)IMTRL, IBARST (J)
         CALL LTSORT (MS, LINKB, IBARST (J), IZ, ADDLNK)
         RETURN
      ELSEIF (IPNTR .EQ. 0) THEN
         IMINUS =  - 1
         CALL LTSORT (MS + MR, LINKM, IMTRL, IMINUS, ADDLNK)
      ENDIF
      JMAT (J) = IMTRL

      RETURN

10000 FORMAT ('   OLD BAR SET NO:', I5, '   TO NEW BAR SET NO:', I5)
10010 FORMAT (' BAR SET:', I5, ' HAS LESS THAN ONE LINE',  / ,
     &   ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
10020 FORMAT (' MATERIAL:', I5, ' FOR BAR SET:', I5,
     &   ' HAS BEEN DESIGNATED',
     &   / , ' AS A REGION  (4 NODE ELEMENT) MATERIAL.',  / ,
     &   ' ELEMENTS WITH 2 AND 4 NODES CANNOT SHARE MATERIAL ID''S',
     &   / , ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
      END
