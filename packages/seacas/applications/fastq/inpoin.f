C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INPOIN (MP, N1, N18, JJ, X, Y, NHOLDP, IHOLDP, IPOINT,
     &   COOR, IPBOUN, LINKP, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE INPOIN = ENTERS A POINT INTO THE DATABASE

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP), LINKP (2, MP)
      DIMENSION IHOLDP (2, MP)

      LOGICAL NOROOM, MERGE, ADDLNK

      NOROOM = .TRUE.
      JHOLD = JJ

C  ZERO OUT THE LINK ARRAY IF NEEDED

      IF (JJ .GT. N18) THEN
         N18 = JJ

C  GET THE CORRECT NODE NUMBER IF MERGING

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

C  INPUT THE POINT DATA

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

10000 FORMAT ('  OLD POINT NO:', I5, '  TO NEW POINT NO:', I5)
      END
