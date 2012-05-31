C $Id: dlpara.f,v 1.1 1990/11/30 11:06:13 gdsjaar Exp $
C $Log: dlpara.f,v $
C Revision 1.1  1990/11/30 11:06:13  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]DLPARA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DLPARA (X1, Y1, X2, Y2, XM, B, BAD)
C***********************************************************************
C
C  SUBROUTINE DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     INREGN = INPUTS REGION CONNECTIVITIES
C
C***********************************************************************
C
C  VARIABLES USED:
C     X1    = X VALUE OF POINT 1
C     X2    = X VALUE OF POINT 2
C     Y1    = Y VALUE OF POINT 1
C     Y2    = Y VALUE OF POINT 2
C     XM    = THE SLOPE OF A STRIGHT LINE BETWEEN POINT 1 AND 2
C     B     = THE Y INTERCEPT OF THE STRAIGHT LINE BETWEEN POINT 1 AND 2
C
C***********************************************************************
C
      LOGICAL BAD
C
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
