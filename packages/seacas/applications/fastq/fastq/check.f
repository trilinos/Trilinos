C $Id: check.f,v 1.1 1990/11/30 11:04:27 gdsjaar Exp $
C $Log: check.f,v $
C Revision 1.1  1990/11/30 11:04:27  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CHECK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHECK (MIN, MAX, I)
C***********************************************************************
C
C  SUBROUTINE CHECK = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     LIST   = LISTS POINTS,  LINES,  AND REGIONS USED IN MESH DEFINITION
C     ERASE  = DELETES POINTS,  LINES,  AND REGIONS FROM THE MESH
C              DEFINITIONS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MIN   = MINIMUM VALUE TO BE TESTED
C     MAX   = MAXIMUM VALUE TO BE TESTED
C     I     = THE ABSOLUTE MAXIMUM VALUE ALLOWED  (THE MINIMUM IS 1)
C
C************************************************************************
C
      IF (MIN .LT. 1)MIN = 1
      IF (MAX .GT. I)MAX = I
      IF (MAX .LT. MIN)MAX = MIN
      RETURN
      END
