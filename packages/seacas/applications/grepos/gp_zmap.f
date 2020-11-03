C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMAP (LIST, LENLST, MAP)
C=======================================================================
      INTEGER LIST(LENLST)
      INTEGER MAP(*)

      DO 10 I=1, LENLST
         LIST(I) = MAP(LIST(I))
   10 CONTINUE

      RETURN
      END
