C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SPAWN (VAXVMS)
C***********************************************************************

C  SUBROUTINE SPAWN = SPAWNS A PROCESS FOR ESCAPE OUT OF FASTQ

C***********************************************************************

C  VARIABLES USED:
C     VAXVMS = .TRUE. IF RUNNING ON A VAXVMS SYSTEM

C***********************************************************************

      LOGICAL VAXVMS

      IF (VAXVMS) THEN
         continue
      ELSE
         CALL MESAGE ('SPAWNING POSSIBLE ONLY ON VAXVMS SYSTEM')
      ENDIF

      END
