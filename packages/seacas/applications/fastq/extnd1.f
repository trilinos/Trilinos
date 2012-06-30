C $Id: extnd1.f,v 1.2 1998/07/14 18:18:50 gdsjaar Exp $
C $Log: extnd1.f,v $
C Revision 1.2  1998/07/14 18:18:50  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:09  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:08  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXTND1.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXTND1 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
C***********************************************************************
C
C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), ANGLE (MXND)
      DIMENSION X(1), Y(1)
C
      CANG = (ANGLE (N2) * .5)
      ANG = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2)) - CANG
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N3) - YN (N2)) **2 +  (XN (N3) - XN (N2)) **2)
      DIST = (DIST1 + DIST2) * .5
      IF (CANG .EQ. 0.) THEN
         ADIST = DIST
      ELSE
         ADIST = DIST / SIN (CANG)
      ENDIF
C
      X(1) = ADIST * COS (ANG) + XN (N2)
      Y(1) = ADIST * SIN (ANG) + YN (N2)
C
      RETURN
C
      END
