C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      LOGICAL FUNCTION SHRUNK (RATIO, NROW)
C***********************************************************************

C  FUNCTION SHRUNK = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ELEMENT
C                    SIZE IS DIMINISHING WITH ROW DEPTH

C***********************************************************************

      DATA TOLER1 /.85/, TOLER2 /.75/, TOLER3 /.6/

      IF ((NROW .GE. 3) .AND. (RATIO .LT. TOLER1)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 2) .AND. (RATIO .LT. TOLER2)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 1) .AND. (RATIO .LT. TOLER3)) THEN
         SHRUNK = .TRUE.
      ELSE
         SHRUNK = .FALSE.
      ENDIF

      RETURN
      END
