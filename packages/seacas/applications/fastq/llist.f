C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LLIST (MS, ML, MAXNL, NS, NL, KNUM, LISTL, ILINE,
     &   ISIDE, NLPS, IFLINE, ILLIST, LCON, ISLIST, LINKS, LINKL, ERR)
C***********************************************************************

C  SUBROUTINE LLIST = PRODUCE LIST OF LINES FOR REGION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     PERIM  = GENERATES PERIMETER OF THE REGION

C***********************************************************************

C     PRODUCE THE LIST OF  (PHYSICAL INDICES OF) LINES FOR  (PHYSICAL)
C     REGION KREG.  THIS LIST IS  (LISTL (I), I=1, NL).
C     *BACKWARDS* SIDES WILL BE REVERSED.
C     IT IS ASSUMED LINES ARE PROPERLY LISTED IN ORDER ON SIDE CARDS.
C     IF THEY ARE NOT,  PERIM WILL DIAGNOSE IT.
C     ERR = .TRUE. IF ERRORS WERE ENCOUNTERED.

C***********************************************************************

      DIMENSION ILINE (ML), LCON (3, ML)
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION ISLIST (NS), LISTL (MAXNL)
      DIMENSION LINKL (2, ML), LINKS (2, MS)

      LOGICAL ERR, ADDLNK

      ERR = .TRUE.
      ADDLNK = .FALSE.
      IS = ISLIST (1)

C  FIRST SIDE

      IF (IS .EQ. 0) THEN
         RETURN
      ELSEIF (IS .LT. 0) THEN
         IFL1 = IABS (IS)
         ILL1 = IABS (IS)
         NEW = 1
         LISTL (NEW) = IFL1
      ELSE
         CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
         I1 = IFLINE (IPNTR)
         I2 = I1 + NLPS (IPNTR) - 1
         IFL1 = ILLIST (I1)
         ILL1 = ILLIST (I2)
         NEW = 0
         DO 100 I = I1, I2
            NEW  =  NEW + 1
            LISTL (NEW) = ILLIST (I)
  100    CONTINUE
      ENDIF
      IF (NS .LE. 1) THEN
         NL = NEW
         ERR = .FALSE.
         RETURN
      ENDIF

C     SECOND SIDE

      IS2 = ISLIST (2)
      IF (IS2 .EQ. 0) THEN
         RETURN
      ELSEIF (IS2 .LT. 0) THEN
         IFL2 = IABS (IS2)
         ILL2 = IABS (IS2)
      ELSE
         CALL LTSORT (MS, LINKS, IS2, IPNTR, ADDLNK)
         I1 = IFLINE (IPNTR)
         I2 = I1 + NLPS (IPNTR) - 1
         IFL2 = ILLIST (I1)
         ILL2 = ILLIST (I2)
      ENDIF

C  DECIDE WHICH END OF SIDE ONE IS THE STARTING POINT

      CALL LTSORT (ML, LINKL, IFL2, IPNTR, ADDLNK)
      K1 = LCON (1, IPNTR)
      K2 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, ILL2, IPNTR, ADDLNK)
      K3 = LCON (1, IPNTR)
      K4 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, IFL1, IPNTR, ADDLNK)
      J1 = LCON (1, IPNTR)
      J2 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, ILL1, IPNTR, ADDLNK)
      J3 = LCON (1, IPNTR)
      J4 = LCON (2, IPNTR)

C  FIRST SIDE IN PROPER ORDER

      IF ( (J3 .EQ. K1) .OR. (J3 .EQ. K2) .OR. (J3 .EQ. K3)
     &   .OR. (J3 .EQ. K4) .OR. (J4 .EQ. K1) .OR. (J4 .EQ. K2)
     &   .OR. (J4 .EQ. K3) .OR. (J4 .EQ. K4)) THEN
         CONTINUE

C  FIRST SIDE NEEDS REVERSED

      ELSEIF ( (J1 .EQ. K1) .OR. (J1 .EQ. K2) .OR. (J1 .EQ. K3)
     &   .OR. (J1 .EQ. K4) .OR. (J2 .EQ. K1) .OR. (J2 .EQ. K2)
     &   .OR. (J2 .EQ. K3) .OR. (J2 .EQ. K4)) THEN
         CALL IREVER (LISTL, NEW)

C  CONNECTIVITY DOES NOT EXIST

      ELSE
         IF (IS2 .GT. 0) THEN
            CALL LTSORT (MS, LINKS, IS2, IPNTR, ADDLNK)
            WRITE ( * , 10000)KNUM, ISIDE (IPNTR)
         ELSE
            CALL LTSORT (ML, LINKL, IABS (IS2), IPNTR, ADDLNK)
            WRITE ( * , 10010)KNUM, ILINE (IPNTR)
         ENDIF
         RETURN
      ENDIF

      NL = NEW
      DO 120 KS = 2, NS
         I = LISTL (NL)
         CALL LTSORT (ML, LINKL, I, IPNTR, ADDLNK)
         J1 = LCON (1, IPNTR)
         J2 = LCON (2, IPNTR)
         IS = ISLIST (KS)

C     ADD NEW LINES TO LIST

         IF (IS .EQ. 0) THEN
            RETURN
         ELSEIF (IS .LT. 0) THEN
            IFL = IABS (IS)
            ILL = IABS (IS)
            NEW = NL + 1
            LISTL (NEW) = IABS (IS)
         ELSE
            CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
            I1 = IFLINE (IPNTR)
            I2 = I1 + NLPS (IPNTR) - 1
            IFL = ILLIST (I1)
            ILL = ILLIST (I2)
            NEW = NL
            DO 110 I = I1, I2
               NEW = NEW + 1
               LISTL (NEW) = ILLIST (I)
  110       CONTINUE
         ENDIF

C     DETERMINE WHETHER THIS SIDE IS BACKWARDS

         CALL LTSORT (ML, LINKL, IFL, IPNTR, ADDLNK)
         K1 = LCON (1, IPNTR)
         K2 = LCON (2, IPNTR)
         CALL LTSORT (ML, LINKL, ILL, IPNTR, ADDLNK)
         K3 = LCON (1, IPNTR)
         K4 = LCON (2, IPNTR)

C  THIS SIDE IS IN PROPER ORDER

         IF ( (J1 .EQ. K1) .OR. (J1 .EQ. K2) .OR. (J2 .EQ. K1)
     &      .OR. (J2 .EQ. K2)) THEN
            CONTINUE

C  THIS SIDE NEEDS REVERSING

         ELSEIF ( (J1 .EQ. K3) .OR. (J1 .EQ. K4) .OR. (J2 .EQ. K3)
     &      .OR. (J2 .EQ. K4)) THEN
            CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
            CALL IREVER (LISTL (NL + 1), NLPS (IPNTR))

C  CONNECTIVITY DOES NOT EXIST

         ELSE
            IF (IS .GT. 0) THEN
               CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
               WRITE ( * , 10000)KNUM, ISIDE (IPNTR)
            ELSE
               CALL LTSORT (ML, LINKL, IABS (IS), IPNTR, ADDLNK)
               WRITE ( * , 10010)KNUM, ILINE (IPNTR)
            ENDIF
            RETURN
         ENDIF

         NL = NEW

  120 CONTINUE

C  SUCCESSFUL LINE LIST GENERATION

      ERR = .FALSE.

      RETURN

10000 FORMAT  (' IN REGION', I5, ', SIDE', I5, ' DOES NOT CONNECT TO',
     &   /, ' THE PREVIOUS SECTION OF THE PERIMETER')
10010 FORMAT  (' IN REGION', I5, ', LINE', I5, ' DOES NOT CONNECT TO',
     &   /, ' THE PREVIOUS SECTION OF THE PERIMETER')

      END
