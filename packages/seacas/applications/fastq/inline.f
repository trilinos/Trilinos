C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INLINE (ML, N2, N19, JJ, LTYP, IP1, IP2, IP3, NN,
     &   FACT, NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON,
     &   ILBOUN, ISBOUN, LINKL, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE INLINE = INPUTS A LINE INTO THE DATABASE

C***********************************************************************

      DIMENSION ILINE (ML), LTYPE (ML), NINT (ML), FACTOR (ML)
      DIMENSION LCON (3, ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML), LINKL (2, ML), IHOLDL (2, ML)

      LOGICAL MERGE, NOROOM, ADDLNK

      NOROOM = .TRUE.
      JHOLD = JJ

C  ADJUST THE COUNTER IF NEEDED

      IF (JJ .GT. N19) THEN
         N19 = JJ

C  GET THE CORRECT LINE NUMBER IF MERGING

      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (ML, LINKL, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDL)NHOLDL = JHOLD
            CALL LTSORT (ML, IHOLDL, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N19 + 1
               N19 = JJ
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (ML, IHOLDL, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF

C  INPUT THE LINE DATA

      N2 = N2 + 1
      J = N2
      IF (J .GT. ML)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (ML, LINKL, JJ, J, ADDLNK)
      ILINE (J) = JJ
      LTYPE (J) = LTYP
      LCON (1, J) = IP1
      LCON (2, J) = IP2
      LCON (3, J) = IP3
      IF (NN .LT. 0) THEN
        NINT(J) = -NN
        WRITE (*, 10010) J
      ELSE
        NINT (J) = NN
      END IF
      FACTOR (J) = FACT
      IF (FACTOR (J) .LE. 0.)FACTOR (J) = 1.
      ILBOUN (J) = 0
      ISBOUN (J) = 0
      NOROOM = .FALSE.
      RETURN

10000 FORMAT ('   OLD LINE NO:', I5, '   TO NEW LINE NO:', I5)
10010 FORMAT ('WARNING: Intervals on line ', I5, ' are negative.',
     &  ' Changed to positive.')
      END
