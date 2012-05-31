C $Id: vinter.f,v 1.1 1990/11/30 11:17:36 gdsjaar Exp $
C $Log: vinter.f,v $
C Revision 1.1  1990/11/30 11:17:36  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]VINTER.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE VINTER (MXND, XN, YN, N1, N2, N3, XOLD, YOLD,
     &   XNEW, YNEW, VCROSS)
C***********************************************************************
C
C  SUBROUTINE VINTER = FINDS WHERE A VECTOR FROM N1 TO N2
C                      INTERSECTS THE VECTOR FROM N3 TO (XOLD, YOLD)
C
C***********************************************************************
C
C  NOTE:  THIS INTERSECTION ROUTINE IS BASED ON AN ALGORITHM GIVEN
C         IN THE BOOK "GEOMETRIC MODELING" BY MICHAEL E. MORTENSON ON
C         PAGES 319 - 320.
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
C
      LOGICAL VCROSS
C
      VCROSS = .FALSE.
C
C  SET UP THE FIRST LINE'S VECTORS (A AND B)
C
      XA = XN (N1)
      YA = YN (N1)
      XB = XN (N2) - XN (N1)
      YB = YN (N2) - YN (N1)
C
C  SET UP THE SECOND LINE'S VECTORS (C AND D)
C
      XC = XN (N3)
      YC = YN (N3)
      XD = XOLD - XN (N3)
      YD = YOLD - YN (N3)
C
C  NOW USE THE VECTORS AND SOLVE FOR W.
C  W IS THE PROPORTION OF THE DISTANCE ALONG THE VECTOR D
C  WHERE THE INTERSECTION OCCURS.  LIKEWISE U IS THE PROPORTIONAL
C  DISTANCE ALONG THE VECTOR B FOR THE INTERSECTION.
C
      DENOM = (YB * XD) - (XB * YD)
C
C  CHECK FOR SPECIAL PARALLEL CASE - THE DENOMINATOR IS EQUAL TO ZERO.
C
      IF (DENOM .NE. 0.) THEN
C
C  GET INTERSECTION LOCATION
C
         W = ( (YC * XB) - (XB * YA) - (XC * YB) + (YB * XA) ) / DENOM
C
C  GET THE U VALUE TO CONFIRM.
C
         IF (XB .NE. 0.) THEN
            U = ( XC + (W * XD) - XA ) / XB
         ELSE
            U = ( YC + (W * YD) - YA ) / YB
         ENDIF
C
C  CALCULATE THE INTERSECTION POINT BASED ON SIMILAR TRIANGLES
C
         XNEW = ( (XA + (XB * U)) + (XC + (XD * W)) ) * .5
         YNEW = ( (YA + (YB * U)) + (YC + (YD * W)) ) * .5
         VCROSS = .TRUE.
      ENDIF
C
      RETURN
C
      END
