C $Id: inline.f,v 1.2 1998/12/08 14:26:04 gdsjaar Exp $
C $Log: inline.f,v $
C Revision 1.2  1998/12/08 14:26:04  gdsjaar
C Detect whether negative line intervals entered. Output warning message
C and fix (make positive).
C
C Upped version to 2.10
C
C Revision 1.1.1.1  1990/11/30 11:09:53  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:09:51  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]INLINE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INLINE (ML, N2, N19, JJ, LTYP, IP1, IP2, IP3, NN,
     &   FACT, NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON,
     &   ILBOUN, ISBOUN, LINKL, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INLINE = INPUTS A LINE INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION ILINE (ML), LTYPE (ML), NINT (ML), FACTOR (ML)
      DIMENSION LCON (3, ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML), LINKL (2, ML), IHOLDL (2, ML)
C
      LOGICAL MERGE, NOROOM, ADDLNK
C
      NOROOM = .TRUE.
      JHOLD = JJ
C
C  ADJUST THE COUNTER IF NEEDED
C
      IF (JJ .GT. N19) THEN
         N19 = JJ
C
C  GET THE CORRECT LINE NUMBER IF MERGING
C
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
C
C  INPUT THE LINE DATA
C
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
C
10000 FORMAT ('   OLD LINE NO:', I5, '   TO NEW LINE NO:', I5)
10010 FORMAT ('WARNING: Intervals on line ', I5, ' are negative.',
     &  ' Changed to positive.')
      END
