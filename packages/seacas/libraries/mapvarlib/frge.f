C Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C======================================================================
      SUBROUTINE FRGE(N,S,L,G)

C**********************************************************************

C Subroutine BS does the back substitution for the solution of the
C local least squares extrapolation technique for element variables
C from their element centroid location to a nodal location.
C The least squares solution is started by a Gauss elimination in
C subroutine FRGE. The process is started in subroutines EXTQ for
C 4-node quads or EXTH for 8-node hexes.

C Called by EXTQ & EXTH

C**********************************************************************

C N     INT   number of equations - 1 + the number of dimensions
C S     REAL  the coefficient matrix
C G     REAL  dummy array
C L     INT   dummy array - placeholder for subscripts
C X     REAL  the solution vector - coefficients of the equation
C SMAX  REAL  dummy variable - used in the solution scheme
C RMAX  REAL  dummy variable - used in the solution scheme
C XMULT REAL  dummy variable - used in the solution scheme
C R     REAL  dummy variable - used in the solution scheme

C**********************************************************************
      DOUBLE PRECISION S(N,N),G(N),SMAX,RMAX,XMULT,R
      INTEGER L(N)

      DO 3 I = 1,N
        L(I) = I
        SMAX = 0.D+00
        DO 2 J = 1,N
          SMAX = MAX(SMAX,DABS(S(I,J)))
    2   CONTINUE
        G(I) = SMAX
    3 CONTINUE

      DO 7 K = 1,N-1
        RMAX = 0.D+00
        JJ = 0
        DO 4 I = K,N
          R = DABS(S(L(I),K)) / G(L(I))
          IF (R .LE. RMAX) GO TO 4
          JJ = I
          RMAX = R
    4   CONTINUE
        if (jj .ne. 0) then
          LK = L(JJ)
          L(JJ) = L(K)
          L(K) = LK
        end if
        DO 6 I = K+1,N
          XMULT = S(L(I),K)/S(L(K),K)
          DO 5 J = K+1,N
            S(L(I),J) = S(L(I),J) - XMULT * S(L(K),J)
    5     CONTINUE
          S(L(I),K) = XMULT
    6   CONTINUE
    7 CONTINUE
      RETURN
      END
