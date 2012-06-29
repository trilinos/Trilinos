C $Id: exdsct.f,v 1.3 1998/07/14 18:18:46 gdsjaar Exp $
C $Log: exdsct.f,v $
C Revision 1.3  1998/07/14 18:18:46  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/03/21 15:44:46  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:07:01  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:07:00  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXDSCT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXDSCT (MXND, XN, YN, LNODES, ANGLE, N1, XNEW, YNEW)
C***********************************************************************
C
C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A DISSECTION NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)
C
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
C
C      XNEW = XN (N0) + XN (N2) - XN (N1)
C      YNEW = YN (N0) + YN (N2) - YN (N1)

      PID2 = 0.5 * ATAN2(0.0, -1.0)

      ANG = ATAN2 (YN (N0) - YN (N1), XN (N0) - XN (N1)) -
     &   (ANGLE (N1) - PID2)
C     &   (2. * ANGLE(N1) / 3.)
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N0) - YN (N1)) **2 +  (XN (N0) - XN (N1)) **2)
      DIST = (DIST1 + DIST2) * .5
      XNEW = DIST * COS (ANG) + XN (N1)
      YNEW = DIST * SIN (ANG) + YN (N1)
C
      RETURN
C
      END
