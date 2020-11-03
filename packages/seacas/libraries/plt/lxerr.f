C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LXERR(MSG,DISP)
      CHARACTER*(*) MSG
      INTEGER DISP
      CHARACTER*80 LOCMSG

      LOCMSG = MSG
      IF (DISP.GE.2) THEN
         CALL LXCLN
      END IF

      CALL SIORPT('LEX',LOCMSG,DISP)
      RETURN

      END
