C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CHECK (MIN, MAX, I)
C***********************************************************************

C  SUBROUTINE CHECK = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     LIST   = LISTS POINTS,  LINES,  AND REGIONS USED IN MESH DEFINITION
C     ERASE  = DELETES POINTS,  LINES,  AND REGIONS FROM THE MESH
C              DEFINITIONS

C***********************************************************************

C  VARIABLES USED:
C     MIN   = MINIMUM VALUE TO BE TESTED
C     MAX   = MAXIMUM VALUE TO BE TESTED
C     I     = THE ABSOLUTE MAXIMUM VALUE ALLOWED  (THE MINIMUM IS 1)

C************************************************************************

      IF (MIN .LT. 1)MIN = 1
      IF (MAX .GT. I)MAX = I
      IF (MAX .LT. MIN)MAX = MIN
      RETURN
      END
