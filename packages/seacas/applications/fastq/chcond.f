C $Id: chcond.f,v 1.1 1990/11/30 11:04:25 gdsjaar Exp $
C $Log: chcond.f,v $
C Revision 1.1  1990/11/30 11:04:25  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CHCOND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHCOND (NPER, NSA, SMANG, INDEX, IFIRST, N1, N2, N3,
     &   N4, CHECK)
C*********************************************************************
C
C  SUBROUTINE CHCOND = THIS SUBROUTINE CHECKS IF THE "NSA" ANGLES
C                       SATISFIES THE CONDITIONS.
C
C*********************************************************************
C
C  VARIABLES  IN: NPER .... NUMBER OF POINTS IN THE REGION
C                 NSA ..... NUMBER OF SOUGHT SMALLEST ANGLES
C                 SMANG ... ARRAY OF SMALLEST ANGLES
C               INDEX ... POINTERS INTO THE ANGLE ARRAY
C           OUT: IFIRST... POINTER TO THE FIRST VERTEX
C                 Mi ...... INTERVALS FOR THE PENTAGON REGION
C                 CHECK ... .EQ. TRUE IF IT SATISFIES THE CONDITIONS
C
C  CALL BY: PICKM5.FOR
C
C WRITTEN BY: HORACIO RECALDE                          DATE:FEB 15, 1988
C
C************************************************************************
C
      PARAMETER (NSANG = 10)
      DIMENSION SMANG(NSA + 1), NAUX(NSANG), INDEX(NSA + 1)
      LOGICAL CHECK
C
      NSA2 = NSA/2
C
C--- SORT THE INDEX ARRAY TO FIND THE 'NSA2' SMALLEST ANGLES
C
      CALL SORTIA (NSA, INDEX, NSA2, NAUX)
      IFIRST = NAUX(1)
      N1 = NAUX(2) - NAUX(1)
      N2 = NAUX(3) - NAUX(2)
      N3 = NAUX(4) - NAUX(3)
      N4 = NAUX(5) - NAUX(4)
      N5 = NPER - N1 - N2 - N3 - N4
C
C--- CHECK COMPATIBILITY EQUATIONS
C
      IF ((N1 + N2 + N3 .LT. N4 + N5 + 2) .OR.
     &   (N2 + N3 + N4 .LT. N5 + N1 + 2) .OR.
     &   (N3 + N4 + N5 .LT. N1 + N2 + 2) .OR.
     &   (N4 + N5 + N1 .LT. N2 + N3 + 2) .OR.
     &   (N5 + N1 + N2 .LT. N3 + N4 + 2)) THEN
         CHECK = .FALSE.
      ELSE
         CHECK = .TRUE.
      ENDIF
C
      RETURN
      END
