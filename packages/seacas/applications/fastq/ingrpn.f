C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INGRPN (MS, MR, N7, N8, N22, JJ, IIN, IFOUND, IREGN,
     &   NSPR, IFSIDE, ISLIST, LINKR, NHOLDR, IHOLDR, IRGFLG, MERGE,
     &   NOROOM)
C***********************************************************************

C  SUBROUTINE INGRPN = INPUTS A REGION INTO THE DATABASE

C***********************************************************************

      DIMENSION IREGN(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION LINKR(2, MR), IHOLDR(2, MR), IRGFLG(MR), IIN(IFOUND)

      LOGICAL NOROOM, MERGE, ADDLNK

      IZ = 0
      NOROOM = .TRUE.
      ADDLNK = .FALSE.

C  ZERO THE LINK ARRAY IF NEEDED

      IF (JJ .GT. N22) THEN
         N22 = JJ

C  SET UP POINTERS FOR MERGING DATA

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

C  ADD THE REGION INTO THE DATABASE

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

      NOROOM = .FALSE.
      RETURN

10000 FORMAT (' OLD GROUP NO:', I5, ' TO NEW GROUP NO:', I5)
10010 FORMAT (' GROUP:', I5, ' HAS LESS THAN ONE REGION', /,
     &   ' THIS GROUP WILL NOT BE INPUT INTO DATABASE')
      END
