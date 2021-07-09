C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INTSCT (X1, Y1, X2, Y2, X3, Y3, X4, Y4, U, W, LCROSS)
C***********************************************************************

C  SUBROUTINE INTSCT = CHECKS TO SEE IF THE LINE FROM N1 TO N2
C                      INTERSECTS THE LINE FROM N3 TO N4

C***********************************************************************

C  NOTE:  THIS INTERSECTION ROUTINE IS BASED ON AN ALGORITHM GIVEN
C         IN THE BOOK "GEOMETRIC MODELING" BY MICHAEL E. MORTENSON ON
C         PAGES 319 - 320.

C***********************************************************************

      LOGICAL LCROSS

      LCROSS = .FALSE.

      if (max(x1, x2) .lt. min(x3, x4)) return
      if (max(y1, y2) .lt. min(y3, y4)) return
      if (max(x3, x4) .lt. min(x1, x2)) return
      if (max(y3, y4) .lt. min(y1, y2)) return

C  SET UP THE FIRST LINE'S VECTORS (A AND B)

      XA = X1
      YA = Y1
      XB = X2 - X1
      YB = Y2 - Y1

C  SET UP THE SECOND LINE'S VECTORS (C AND D)

      XC = X3
      YC = Y3
      XD = X4 - X3
      YD = Y4 - Y3

C  NOW USE THE VECTORS AND SOLVE FOR W.
C  W IS THE PROPORTION OF THE DISTANCE ALONG THE VECTOR D
C  WHERE THE INTERSECTION OCCURS.  LIKEWISE U IS THE PROPORTIONAL
C  DISTANCE ALONG THE VECTOR B FOR THE INTERSECTION.   IF THERE IS
C  AN INTERSECTION, BOTH U AND W MUST BE BETWEEN 0 AND 1.

      DENOM = (YB * XD) - (XB * YD)

C  CHECK FOR SPECIAL PARALLEL CASE - THE DENOMINATOR IS EQUAL TO ZERO.

      IF (DENOM .NE. 0.) THEN

C  CHECK FOR INTERSECTION

         W = ( (YC * XB) - (XB * YA) - (XC * YB) + (YB * XA) ) / DENOM
         IF ( (W .LT. 1.) .AND. (W .GT. 0.) ) THEN

C  W INDICATES AN INTERSECTION HAS OCCURRED.
C  GET THE U VALUE AND CONFIRM.

            IF (XB .NE. 0.) THEN
               U = ( XC + (W * XD) - XA ) / XB
            ELSE
               U = ( YC + (W * YD) - YA ) / YB
            ENDIF
            IF ( (U .LT. 1.) .AND. (U .GT. 0.) ) THEN
               LCROSS = .TRUE.
            ENDIF
         ENDIF
      ENDIF

      RETURN

      END
