C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INCMAP (IDMAP, NUMTYPE, IDINC)
C=======================================================================

      INTEGER       IDMAP(*)

C ... This routine increments all IDS by IDINC.

      IF (NUMTYPE .LE. 0) RETURN
      DO 20 ID = 1, NUMTYPE
         IDMAP(ID) = IDMAP(ID) + IDINC
   20 CONTINUE

      RETURN
      END
