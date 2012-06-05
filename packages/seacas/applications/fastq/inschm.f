C $Id: inschm.f,v 1.1 1990/11/30 11:10:12 gdsjaar Exp $
C $Log: inschm.f,v $
C Revision 1.1  1990/11/30 11:10:12  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INSCHM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INSCHM (MR, MSC, N8, N19, JJ, DUMMY, ISCHM, SCHEME,
     &   LINKSC, DEFSCH, NOROOM, DOLINK)
C***********************************************************************
C
C  SUBROUTINE INSCHM = INPUTS A SCHEME INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION ISCHM (MSC), SCHEME (MSC), LINKSC (2, MR)
C
      CHARACTER * 72 SCHEME, DEFSCH, DUMMY
C
      LOGICAL NOROOM, DOLINK, ADDLNK
C
      NOROOM = .TRUE.
      ADDLNK = .TRUE.
C
C  ENTER THE DEFAULT SCHEME IF THE REGION NUMBER IS ZERO
C
      IF (JJ .EQ. 0) THEN
         DEFSCH = DUMMY
C
C  ENTER THE SCHEME
C
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
C
      END
