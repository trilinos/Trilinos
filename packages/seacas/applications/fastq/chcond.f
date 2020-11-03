C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CHCOND (NPER, NSA, SMANG, INDEX, IFIRST, N1, N2, N3,
     &   N4, CHECK)
C*********************************************************************

C  SUBROUTINE CHCOND = THIS SUBROUTINE CHECKS IF THE "NSA" ANGLES
C                       SATISFIES THE CONDITIONS.

C*********************************************************************

C  VARIABLES  IN: NPER .... NUMBER OF POINTS IN THE REGION
C                 NSA ..... NUMBER OF SOUGHT SMALLEST ANGLES
C                 SMANG ... ARRAY OF SMALLEST ANGLES
C               INDEX ... POINTERS INTO THE ANGLE ARRAY
C           OUT: IFIRST... POINTER TO THE FIRST VERTEX
C                 Mi ...... INTERVALS FOR THE PENTAGON REGION
C                 CHECK ... .EQ. TRUE IF IT SATISFIES THE CONDITIONS

C  CALL BY: PICKM5.FOR

C WRITTEN BY: HORACIO RECALDE                          DATE:FEB 15, 1988

C************************************************************************

      PARAMETER (NSANG = 10)
      DIMENSION SMANG(NSA + 1), NAUX(NSANG), INDEX(NSA + 1)
      LOGICAL CHECK

      NSA2 = NSA/2

C--- SORT THE INDEX ARRAY TO FIND THE 'NSA2' SMALLEST ANGLES

      CALL SORTIA (NSA, INDEX, NSA2, NAUX)
      IFIRST = NAUX(1)
      N1 = NAUX(2) - NAUX(1)
      N2 = NAUX(3) - NAUX(2)
      N3 = NAUX(4) - NAUX(3)
      N4 = NAUX(5) - NAUX(4)
      N5 = NPER - N1 - N2 - N3 - N4

C--- CHECK COMPATIBILITY EQUATIONS

      IF ((N1 + N2 + N3 .LT. N4 + N5 + 2) .OR.
     &   (N2 + N3 + N4 .LT. N5 + N1 + 2) .OR.
     &   (N3 + N4 + N5 .LT. N1 + N2 + 2) .OR.
     &   (N4 + N5 + N1 .LT. N2 + N3 + 2) .OR.
     &   (N5 + N1 + N2 .LT. N3 + N4 + 2)) THEN
         CHECK = .FALSE.
      ELSE
         CHECK = .TRUE.
      ENDIF

      RETURN
      END
