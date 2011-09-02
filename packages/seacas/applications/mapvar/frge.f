C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C======================================================================
      SUBROUTINE FRGE(N,S,L,G)
C
C**********************************************************************
C
C Subroutine BS does the back substitution for the soultion of the 
C local least squares extrapolation technique for element variables
C from their element centroid location to a nodal location.
C The least squares solution is started by a Gauss elimination in
C subroutine FRGE. The process is started in subroutines EXTQ for
C 4-node quads or EXTH for 8-node hexes.
C
C Called by EXTQ & EXTH
C
C**********************************************************************
C
C N     INT   number of equations - 1 + the number of dimensions
C S     REAL  the coefficient matrix
C G     REAL  dummy array
C L     INT   dummy array - placeholder for subscripts
C X     REAL  the solution vector - coefficients of the equation
C SMAX  REAL  dummy variable - used in the solution scheme
C RMAX  REAL  dummy variable - used in the solution scheme
C XMULT REAL  dummy variable - used in the solution scheme
C R     REAL  dummy variable - used in the solution scheme
C
C**********************************************************************
      DOUBLE PRECISION S(N,N),G(N),SMAX,RMAX,XMULT,R
      INTEGER L(N)
C
C
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
        DO 4 I = K,N
          R = DABS(S(L(I),K)) / G(L(I))
          IF (R .LE. RMAX) GO TO 4
          J = I
          RMAX = R
    4   CONTINUE
        LK = L(J)
        L(J) = L(K)
        L(K) = LK
        DO 6 I = K+1,N
          XMULT = S(L(I),K)/S(LK,K)
          DO 5 J = K+1,N
            S(L(I),J) = S(L(I),J) - XMULT * S(LK,J)
    5     CONTINUE
          S(L(I),K) = XMULT
    6   CONTINUE
    7 CONTINUE
      RETURN
      END
