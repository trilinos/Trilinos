C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DMESS (DEV1, TEXT)
C***********************************************************************

C  SUBROUTINE DMESS = PRINTS A ONE LINE MESSAGE AT THE BOTTOM OF THE
C                       SCREEN

C***********************************************************************

      CHARACTER*(*) TEXT, DEV1*3

      CALL MESAGE (TEXT)
      RETURN

      END
