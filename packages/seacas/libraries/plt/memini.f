C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MEMINI(MEMRY,LENGTH)
      INTEGER MEMRY(*)
      INTEGER LENGTH

      DO 2630 I = 1,LENGTH
         MEMRY(1) = LENGTH
 2630 CONTINUE
      RETURN

      END
