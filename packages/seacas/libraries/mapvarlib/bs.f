C Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C======================================================================
      SUBROUTINE BS(N,S,F,L,X)

C**********************************************************************

C Subroutine BS does the back substitution for the solution of the
C local least squares extrapolation technique for element variables
C from their element centroid location to a nodal location.
C The least squares solution is started by a Gauss elimination in
C subroutine FRGE. The process is started in subroutines EXTQ for
C 4-node quads or EXTH for 8-node hexes.

C Called by EXTH, EXTQ

C**********************************************************************

C N    INT   number of equations - 1 + the number of dimensions
C S    REAL  the coefficient matrix - after forward gauss elimination
C F    REAL  the load vector
C L    INT   dummy array - placeholder for subscripts
C X    REAL  the solution vector - coefficients of the equation
C SUM  REAL  dummy variable - used in the solution scheme

C**********************************************************************

C subroutine written in double precision

      DOUBLE PRECISION S(N,N),F(N),X(N),SUM
      INTEGER L(N)

      DO 3 K = 1,N-1
        DO 2 I = K+1,N
          F(L(I)) = F(L(I)) - S(L(I),K) * F(L(K))
    2   CONTINUE
    3 CONTINUE
      X(N) = F(L(N)) / S(L(N),N)
      DO 5 I = N-1,1,-1
        SUM = F(L(I))
        DO 4 J = I+1,N
          SUM = SUM - S(L(I),J) * X(J)
    4   CONTINUE
        X(I) = SUM / S(L(I),I)
    5 CONTINUE
      RETURN
      END
