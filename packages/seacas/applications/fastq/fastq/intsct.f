C $Id: intsct.f,v 1.2 1992/02/04 15:50:25 gdsjaar Exp $
C $Log: intsct.f,v $
C Revision 1.2  1992/02/04 15:50:25  gdsjaar
C Added bouding box check in intsct, reduce time 5 percent
C
c Revision 1.1.1.1  1990/11/30  11:10:21  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:10:20  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]INTSCT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INTSCT (X1, Y1, X2, Y2, X3, Y3, X4, Y4, U, W, LCROSS)
C***********************************************************************
C
C  SUBROUTINE INTSCT = CHECKS TO SEE IF THE LINE FROM N1 TO N2
C                      INTERSECTS THE LINE FROM N3 TO N4
C
C***********************************************************************
C
C  NOTE:  THIS INTERSECTION ROUTINE IS BASED ON AN ALGORITHM GIVEN
C         IN THE BOOK "GEOMETRIC MODELING" BY MICHAEL E. MORTENSON ON
C         PAGES 319 - 320.
C
C***********************************************************************
C
      LOGICAL LCROSS
C
      LCROSS = .FALSE.

      if (max(x1, x2) .lt. min(x3, x4)) return
      if (max(y1, y2) .lt. min(y3, y4)) return
      if (max(x3, x4) .lt. min(x1, x2)) return
      if (max(y3, y4) .lt. min(y1, y2)) return
C
C  SET UP THE FIRST LINE'S VECTORS (A AND B)
C
      XA = X1
      YA = Y1
      XB = X2 - X1
      YB = Y2 - Y1
C
C  SET UP THE SECOND LINE'S VECTORS (C AND D)
C
      XC = X3
      YC = Y3
      XD = X4 - X3
      YD = Y4 - Y3
C
C  NOW USE THE VECTORS AND SOLVE FOR W.
C  W IS THE PROPORTION OF THE DISTANCE ALONG THE VECTOR D
C  WHERE THE INTERSECTION OCCURS.  LIKEWISE U IS THE PROPORTIONAL
C  DISTANCE ALONG THE VECTOR B FOR THE INTERSECTION.   IF THERE IS
C  AN INTERSECTION, BOTH U AND W MUST BE BETWEEN 0 AND 1.
C
      DENOM = (YB * XD) - (XB * YD)
C
C  CHECK FOR SPECIAL PARALLEL CASE - THE DENOMINATOR IS EQUAL TO ZERO.
C
      IF (DENOM .NE. 0.) THEN
C
C  CHECK FOR INTERSECTION
C
         W = ( (YC * XB) - (XB * YA) - (XC * YB) + (YB * XA) ) / DENOM
         IF ( (W .LT. 1.) .AND. (W .GT. 0.) ) THEN
C
C  W INDICATES AN INTERSECTION HAS OCCURRED.
C  GET THE U VALUE AND CONFIRM.
C
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
C
      RETURN
C
      END
