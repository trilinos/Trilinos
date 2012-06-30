C $Id: inpoin.f,v 1.1 1990/11/30 11:09:57 gdsjaar Exp $
C $Log: inpoin.f,v $
C Revision 1.1  1990/11/30 11:09:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INPOIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INPOIN (MP, N1, N18, JJ, X, Y, NHOLDP, IHOLDP, IPOINT,
     &   COOR, IPBOUN, LINKP, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INPOIN = ENTERS A POINT INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP), LINKP (2, MP)
      DIMENSION IHOLDP (2, MP)
C
      LOGICAL NOROOM, MERGE, ADDLNK
C
      NOROOM = .TRUE.
      JHOLD = JJ
C
C  ZERO OUT THE LINK ARRAY IF NEEDED
C
      IF (JJ .GT. N18) THEN
         N18 = JJ
C
C  GET THE CORRECT NODE NUMBER IF MERGING
C
      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (MP, LINKP, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDP)NHOLDP = JHOLD
            CALL LTSORT (MP, IHOLDP, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N18 + 1
               N18 = N18 + 1
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MP, IHOLDP, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  INPUT THE POINT DATA
C
      N1 = N1 + 1
      J = N1
      IF (J .GT. MP)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MP, LINKP, JJ, J, ADDLNK)
      IPOINT (J) = JJ
      COOR (1, J) = X
      COOR (2, J) = Y
      IPBOUN (J) = 0
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT ('  OLD POINT NO:', I5, '  TO NEW POINT NO:', I5)
      END
