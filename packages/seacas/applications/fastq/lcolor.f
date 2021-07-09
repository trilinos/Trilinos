C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LCOLOR (COLOR)
C***********************************************************************

C  SUBROUTINE LCOLOR = SETS THE LINE COLOR

C***********************************************************************

      CHARACTER*5 COLOR

      IF (COLOR .EQ. 'WHITE') THEN
         CALL PLTSTD(1, 7.)
      ELSEIF (COLOR .EQ. 'BLACK') THEN
         CALL PLTSTD(1, 0.)
      ELSEIF (COLOR .EQ. 'YELOW') THEN
         CALL PLTSTD(1, 3.)
      ELSEIF (COLOR .EQ. 'PINK ') THEN
         CALL PLTSTD(1, 5.)
      ELSEIF (COLOR .EQ. 'RED  ') THEN
         CALL PLTSTD(1, 1.)
      ELSEIF (COLOR .EQ. 'BLUE ') THEN
         CALL PLTSTD(1, 4.)
      ENDIF
      RETURN

      END
