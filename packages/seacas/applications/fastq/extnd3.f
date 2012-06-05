C $Id: extnd3.f,v 1.2 1998/07/14 18:18:51 gdsjaar Exp $
C $Log: extnd3.f,v $
C Revision 1.2  1998/07/14 18:18:51  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:12  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:11  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXTND3.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXTND3 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
C***********************************************************************
C
C  SUBROUTINE EXTND3 = CALCULATES TWO POSITIONS AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE AND ONE AT 1/3 ANGLE INTERVALS
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
      DIMENSION X(3), Y(3)
C
      ANG = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2))
      CANG = ANGLE / 3.
      ANG1 = ANG - CANG
      ANG2 = ANG - (ANGLE * .5)
      ANG3 = ANG - 2. * CANG
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N3) - YN (N2)) **2 +  (XN (N3) - XN (N2)) **2)
      DIST = (DIST1 + DIST2) * .5
      IF (CANG. EQ. 0.) THEN
         ADIST = DIST
      ELSE
         ADIST = DIST / SIN (CANG)
      ENDIF
C
      X(1) = ADIST * COS (ANG1) + XN (N2)
      Y(1) = ADIST * SIN (ANG1) + YN (N2)
      X(2) = 1.4142 * ADIST * COS (ANG2) + XN (N2)
      Y(2) = 1.4142 * ADIST * SIN (ANG2) + YN (N2)
      X(3) = ADIST * COS (ANG3) + XN (N2)
      Y(3) = ADIST * SIN (ANG3) + YN (N2)
C
      RETURN
C
      END
