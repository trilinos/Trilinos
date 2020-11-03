C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DLPARA (X1, Y1, X2, Y2, XM, B, BAD)
C***********************************************************************

C  SUBROUTINE DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     INREGN = INPUTS REGION CONNECTIVITIES

C***********************************************************************

C  VARIABLES USED:
C     X1    = X VALUE OF POINT 1
C     X2    = X VALUE OF POINT 2
C     Y1    = Y VALUE OF POINT 1
C     Y2    = Y VALUE OF POINT 2
C     XM    = THE SLOPE OF A STRIGHT LINE BETWEEN POINT 1 AND 2
C     B     = THE Y INTERCEPT OF THE STRAIGHT LINE BETWEEN POINT 1 AND 2

C***********************************************************************

      LOGICAL BAD

      IF (ABS (X2 - X1) .LT. 0.000001) THEN
         BAD = .TRUE.
         B = X1
      ELSE
         BAD = .FALSE.
         XM =  (Y2 - Y1) / (X2 - X1)
         B = Y1- (X1 * XM)
      ENDIF
      RETURN
      END
