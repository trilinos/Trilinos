C $Id: inregn.f,v 1.2 2004/01/21 05:18:40 gdsjaar Exp $
C $Log: inregn.f,v $
C Revision 1.2  2004/01/21 05:18:40  gdsjaar
C Initialized several variables identified by valgrind.
C
C Revision 1.1.1.1  1990/11/30 11:10:08  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:10:06  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]INREGN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INREGN (MS, MR, N7, N8, N22, N23, JJ, JMTRL, IIN,
     &   IFOUND, IREGN, IMAT, NSPR, IFSIDE, ISLIST, LINKR, LINKM,
     &   NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INREGN = INPUTS A REGION INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION LINKR(2, MR), LINKM(2, (MS + MR))
      DIMENSION IHOLDR(2, MR), IHOLDM(2, (MS + MR)), IRGFLG(MR)
      DIMENSION IIN(IFOUND)
C
      LOGICAL NOROOM, MERGE, ADDLNK
C
      IPNTR = 0
      IZ = 0
      NOROOM = .TRUE.
      ADDLNK = .FALSE.
      IMTRL = ABS(JMTRL)
C
C  ZERO THE LINK ARRAY IF NEEDED
C
      IF (JJ .GT. N22) THEN
         N22 = JJ
C
C  SET UP POINTERS FOR MERGING DATA
C
      ELSE IF (MERGE) THEN
         JHOLD = JJ
         CALL LTSORT (MR, LINKR, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDR)NHOLDR = JHOLD
            CALL LTSORT (MR, IHOLDR, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N22 + 1
               N22 = JJ
               WRITE(*, 10000) JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MR, IHOLDR, JHOLD, JJ, ADDLNK)
            END IF
         END IF
      END IF
C
C  ADD THE REGION INTO THE DATABASE
C
      N7 = N7 + 1
      J = N7
      IF (J .GT. MR) RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MR, LINKR, JJ, J, ADDLNK)
      IREGN(J) = JJ
      IFSIDE(J) = N8 + 1
      IRGFLG(J) = -1
      DO 100 I = 1, IFOUND
         JJ = IIN(I)
         IF (JJ .EQ. 0) GO TO 110
         N8 = N8 + 1
         IF (N8 .GT. MR*4) RETURN
         ISLIST(N8) = JJ
  100 CONTINUE
  110 CONTINUE
      NSPR(J) = N8 - IFSIDE(J) + 1
      IF (NSPR(J) .LT. 1) THEN
         WRITE(*, 10010) J
         CALL LTSORT (MR, LINKR, IREGN(J), IZ, ADDLNK)
      END IF
C
C  LINK UP THE MATERIAL
C
C  ZERO THE LINK ARRAY IF NEEDED
C
      IF (IMTRL .GT. N23) THEN
         N23 = IMTRL
C
C  SET UP POINTERS FOR MERGING DATA
C
      ELSE IF (MERGE) THEN
         JHOLD = IMTRL
         ADDLNK = .FALSE.
         CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
         IF (IPNTR .NE. 0) THEN
            IF (JHOLD .GT. NHOLDM)NHOLDM = JHOLD
            CALL LTSORT (MS + MR, IHOLDM, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               IMTRL = IPNTR
            ELSE
               IMTRL = N23 + 1
               N23 = N23 + 1
               WRITE(*, 10020) JHOLD, IMTRL
               ADDLNK = .TRUE.
               CALL LTSORT (MS + MR, IHOLDM, JHOLD, IMTRL, ADDLNK)
            END IF
         END IF
      END IF
C
C  ADD THE MATERIAL INTO THE DATABASE
C
      NOROOM = .FALSE.
      ADDLNK = .FALSE.
      CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
      IF (IPNTR .LT. 0) THEN
         CALL MESAGE(' ')
         WRITE(*, 10030) IMTRL, IREGN(J)
         ADDLNK = .TRUE.
         CALL LTSORT (MR, LINKR, IREGN(J), IZ, ADDLNK)
         RETURN
      ELSE IF (IPNTR .EQ. 0) THEN
         ADDLNK = .TRUE.
         IONE = 1
         CALL LTSORT (MS + MR, LINKM, IMTRL, IONE, ADDLNK)
      END IF
      IMAT(J) = JMTRL
C
      RETURN
C
10000 FORMAT(' OLD REGION NO:', I5, ' TO NEW REGION NO:', I5)
10010 FORMAT(' REGION:', I5, ' HAS LESS THAN ONE SIDE', /,
     &   ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
10020 FORMAT(' OLD MATERIAL NO:', I5, ' TO NEW MATERIAL NO:', I5)
10030 FORMAT(' MATERIAL:', I5, ' FOR REGION:', I5, ' HAS BEEN '//
     &   'DESIGNATED', /,' AS A BAR SET (2 NODE ELEMENT) MATERIAL.', /,
     &   ' ELEMENTS WITH 2 AND 4 NODES CANNOT SHARE MATERIAL ID''S',/,
     &   ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
      END
