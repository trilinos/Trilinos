C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ERASEC (OLDCUR)
C***********************************************************************

C  SUBROUTINE ERASEC = DEACTIVATES THE CROSSHAIRS

C***********************************************************************

      LOGICAL OLDCUR

      IF (OLDCUR) THEN
         WRITE (*,*) CHAR(27)//'G0'
         OLDCUR = .FALSE.
      ENDIF
      RETURN

      END
