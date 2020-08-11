C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INSCHM (MR, MSC, N8, N19, JJ, DUMMY, ISCHM, SCHEME,
     &   LINKSC, DEFSCH, NOROOM, DOLINK)
C***********************************************************************

C  SUBROUTINE INSCHM = INPUTS A SCHEME INTO THE DATABASE

C***********************************************************************

      DIMENSION ISCHM (MSC), SCHEME (MSC), LINKSC (2, MR)

      CHARACTER * 72 SCHEME, DEFSCH, DUMMY

      LOGICAL NOROOM, DOLINK, ADDLNK

      NOROOM = .TRUE.
      ADDLNK = .TRUE.

C  ENTER THE DEFAULT SCHEME IF THE REGION NUMBER IS ZERO

      IF (JJ .EQ. 0) THEN
         DEFSCH = DUMMY

C  ENTER THE SCHEME

      ELSE
         IF ( (DOLINK) .AND. (JJ .GT. N19))N19 = JJ
         N8 = N8 + 1
         J = N8
         IF (N8 .GT. MSC)RETURN
         ISCHM (J) = JJ
         IF (DOLINK)CALL LTSORT (MR, LINKSC, JJ, J, ADDLNK)
         SCHEME (J) = DUMMY
      ENDIF
      NOROOM = .FALSE.
      RETURN

      END
