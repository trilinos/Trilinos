C $Id: excorn.f,v 1.2 1991/03/21 15:44:44 gdsjaar Exp $
C $Log: excorn.f,v $
C Revision 1.2  1991/03/21 15:44:44  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:06:58  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:57  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXCORN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXCORN (MXND, XN, YN, LNODES, ANGLE, N0, N1, N2, XNEW,
     &   YNEW)
C***********************************************************************
C
C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)
C
      LOGICAL SIDEP
C
C      XNEW = XN (N0) + XN (N2) - XN (N1)
C      YNEW = YN (N0) + YN (N2) - YN (N1)

      PID2 = 0.5 * ATAN2(0.0, -1.0)
C
C      ANG2 = ATAN2 (YN (N1)-YN (N0), XN (N1)-XN (N0))+PID2
      BANG1 = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2))
      BANG2 = ATAN2 (YN (N1) - YN (N0), XN (N1) - XN (N0))
      IF (SIDEP (ANGLE (N2)))THEN
         ANG1 = BANG1 - (ANGLE (N2) * .5)
      ELSE
         ANG1 = BANG1 - PID2
      ENDIF
      IF (SIDEP (ANGLE (N0)))THEN
         ANG2 = BANG2 + (ANGLE (N0) * .5)
      ELSE
         ANG2 = BANG2 + PID2
      ENDIF
      DIST1 = SQRT ((YN (N2)-YN (N1)) **2 +  (XN (N2)-XN (N1)) **2)
      DIST2 = SQRT ((YN (N0)-YN (N1)) **2 +  (XN (N0)-XN (N1)) **2)
      DIST = (DIST1 + DIST2) * .5
      XNEW = ( (DIST * COS (ANG1) + XN (N2)) +
     &   (DIST * COS (ANG2) + XN (N0)) ) * .5
      YNEW = ( (DIST * SIN (ANG1) + YN (N2)) +
     &   (DIST * SIN (ANG2) + YN (N0)) ) * .5
C      XNEW = DIST * COS (ANG1) + XN (N2)
C      YNEW = DIST * SIN (ANG1) + YN (N2)
C
      RETURN
C
      END
