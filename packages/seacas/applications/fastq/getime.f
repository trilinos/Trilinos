C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETIME (TIME)
C***********************************************************************

C  SUBROUTINE GETIME = GETS THE CPU TIME USED BY THE CURRENT PROCESS

C***********************************************************************

      CALL EXCPUS (TIME)
      RETURN
      END
