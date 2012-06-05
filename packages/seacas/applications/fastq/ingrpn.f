C $Id: ingrpn.f,v 1.1 1990/11/30 11:09:41 gdsjaar Exp $
C $Log: ingrpn.f,v $
C Revision 1.1  1990/11/30 11:09:41  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INGRPN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INGRPN (MS, MR, N7, N8, N22, JJ, IIN, IFOUND, IREGN,
     &   NSPR, IFSIDE, ISLIST, LINKR, NHOLDR, IHOLDR, IRGFLG, MERGE,
     &   NOROOM)
C***********************************************************************
C
C  SUBROUTINE INGRPN = INPUTS A REGION INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IREGN(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION LINKR(2, MR), IHOLDR(2, MR), IRGFLG(MR), IIN(IFOUND)
C
      LOGICAL NOROOM, MERGE, ADDLNK
C
      IZ = 0
      NOROOM = .TRUE.
      ADDLNK = .FALSE.
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
            IF (JHOLD .GT. NHOLDR) NHOLDR = JHOLD
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
      IRGFLG(J) = 1
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
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT (' OLD GROUP NO:', I5, ' TO NEW GROUP NO:', I5)
10010 FORMAT (' GROUP:', I5, ' HAS LESS THAN ONE REGION', /,
     &   ' THIS GROUP WILL NOT BE INPUT INTO DATABASE')
      END
