C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE VINTER (MXND, XN, YN, N1, N2, N3, XOLD, YOLD,
     &   XNEW, YNEW, VCROSS)
C***********************************************************************

C  SUBROUTINE VINTER = FINDS WHERE A VECTOR FROM N1 TO N2
C                      INTERSECTS THE VECTOR FROM N3 TO (XOLD, YOLD)

C***********************************************************************

C  NOTE:  THIS INTERSECTION ROUTINE IS BASED ON AN ALGORITHM GIVEN
C         IN THE BOOK "GEOMETRIC MODELING" BY MICHAEL E. MORTENSON ON
C         PAGES 319 - 320.

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)

      LOGICAL VCROSS

      VCROSS = .FALSE.

C  SET UP THE FIRST LINE'S VECTORS (A AND B)

      XA = XN (N1)
      YA = YN (N1)
      XB = XN (N2) - XN (N1)
      YB = YN (N2) - YN (N1)

C  SET UP THE SECOND LINE'S VECTORS (C AND D)

      XC = XN (N3)
      YC = YN (N3)
      XD = XOLD - XN (N3)
      YD = YOLD - YN (N3)

C  NOW USE THE VECTORS AND SOLVE FOR W.
C  W IS THE PROPORTION OF THE DISTANCE ALONG THE VECTOR D
C  WHERE THE INTERSECTION OCCURS.  LIKEWISE U IS THE PROPORTIONAL
C  DISTANCE ALONG THE VECTOR B FOR THE INTERSECTION.

      DENOM = (YB * XD) - (XB * YD)

C  CHECK FOR SPECIAL PARALLEL CASE - THE DENOMINATOR IS EQUAL TO ZERO.

      IF (DENOM .NE. 0.) THEN

C  GET INTERSECTION LOCATION

         W = ( (YC * XB) - (XB * YA) - (XC * YB) + (YB * XA) ) / DENOM

C  GET THE U VALUE TO CONFIRM.

         IF (XB .NE. 0.) THEN
            U = ( XC + (W * XD) - XA ) / XB
         ELSE
            U = ( YC + (W * YD) - YA ) / YB
         ENDIF

C  CALCULATE THE INTERSECTION POINT BASED ON SIMILAR TRIANGLES

         XNEW = ( (XA + (XB * U)) + (XC + (XD * W)) ) * .5
         YNEW = ( (YA + (YB * U)) + (YC + (YD * W)) ) * .5
         VCROSS = .TRUE.
      ENDIF

      RETURN

      END
