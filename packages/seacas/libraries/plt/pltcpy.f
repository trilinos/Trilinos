C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCPY
      INTEGER SUPPO

      CALL PLTFLU
      CALL VDIQES(100,SUPPO)
      IF (SUPPO.NE.0) THEN
         CALL VDESCP(100,0)
      END IF

      RETURN

      END
