C $Id: exside.f,v 1.1 1990/11/30 11:07:05 gdsjaar Exp $
C $Log: exside.f,v $
C Revision 1.1  1990/11/30 11:07:05  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]EXSIDE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXSIDE (MXND, XN, YN, LNODES, ANGLE, N1, XNEW, YNEW)
C***********************************************************************
C
C  SUBROUTINE EXSIDE = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A SIDE NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)
C
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
C
      ANG = ATAN2 (YN (N2) - YN (N1),  XN (N2) - XN (N1)) +
     &   (ANGLE (N1) * .5)
C
      DIST1 = SQRT ( (YN (N2) - YN (N1) ) **2 +
     &   ( XN (N2) - XN (N1) ) **2)
      DIST2 = SQRT ( (YN (N0) - YN (N1) ) **2 +
     &   ( XN (N0) - XN (N1) ) **2)
C
      DIST =  (DIST1 + DIST2) * .5
      XNEW =  (DIST * COS (ANG) ) + XN (N1)
      YNEW =  (DIST * SIN (ANG) ) + YN (N1)
C
      RETURN
C
      END
