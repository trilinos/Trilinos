C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPUMVU(A,B,L)
      IMPLICIT INTEGER (A-Z)
      DIMENSION A(*),B(*)

      DO 2240 I = 1,L
         B(I) = A(I)
 2240 CONTINUE
      RETURN

      END
