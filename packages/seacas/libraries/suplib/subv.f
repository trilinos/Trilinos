C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE SUBV( N,A,B,C )

C***********************************************************************

C     DESCRIPTION: This routine subtracts two vectors

C     FORMAL PARAMETERS:
C        N        INTEGER   Number of entries in A, B
C        A        REAL      First vector
C        B        REAL      Vector to be subtracted
C        C        REAL      Vector with the result

C***********************************************************************

      DIMENSION A(N),B(N),C(N)

      DO 100 I = 1,N
        C(I) = A(I) - B(I)
  100 CONTINUE

      RETURN
      END
