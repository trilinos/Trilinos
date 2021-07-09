C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RGDSTR (NPNODE, NPELEM, KKK, NNXK, XN, YN, NXK)
C************************************************************************

C  SUBROUTINE RGDSTR = CALCULATES A REGION DISTORTION MEASURE

C***********************************************************************

      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)

C  CALCULATE THE ELEMENT DISTORTION

      N1 = NXK (1,1)
      N2 = NXK (2,1)
      N3 = NXK (3,1)
      N4 = NXK (4,1)
      CALL DSTORT (XN (N1), XN (N2), XN (N3), XN (N4),
     &   YN (N1), YN (N2), YN (N3), YN (N4), VALUE)
      VMIN = VALUE
      VMAX = VALUE
      SUM = VALUE
      KMIN = 1
      KMAX = 1
      DO 100 I = 1, KKK
         N1 = NXK (1,I)
         N2 = NXK (2,I)
         N3 = NXK (3,I)
         N4 = NXK (4,I)
         CALL DSTORT (XN (N1), XN (N2), XN (N3), XN (N4),
     &      YN (N1), YN (N2), YN (N3), YN (N4), VALUE)
         IF (VMIN .GT. VALUE) THEN
            VMIN = VALUE
            KMIN = I
         ELSE IF (VMAX .LT. VALUE) THEN
            VMAX = VALUE
            KMAX = I
         ENDIF
         SUM = SUM + VALUE
  100 CONTINUE

C  PRINT OUT THE RESULTS

      SUM = SUM / DBLE(KKK)
      WRITE (*, 10000) VMIN, KMIN, VMAX, KMAX, SUM

      RETURN

10000 FORMAT ('  THE MINIMUM DISTORTION IS: ',G14.7,' IN ELEMENT: ',I10,
     &   /,
     &   '  THE MAXIMUM DISTORTION IS: ',G14.7,' IN ELEMENT: ',I10, /,
     &   '  THE AVERAGE DISTORTION IS: ',G14.7)

      END
