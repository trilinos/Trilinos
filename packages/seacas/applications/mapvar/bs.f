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
      SUBROUTINE BS(N,S,F,L,X)
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
C Called by EXTH, EXTQ
C
C**********************************************************************
C
C N    INT   number of equations - 1 + the number of dimensions
C S    REAL  the coefficient matrix - after forward gauss elimination
C F    REAL  the load vector
C L    INT   dummy array - placeholder for subscripts
C X    REAL  the solution vector - coefficients of the equation
C SUM  REAL  dummy variable - used in the solution scheme
C
C**********************************************************************
C
C subroutine written in double precision
C
      DOUBLE PRECISION S(N,N),F(N),X(N),SUM
      INTEGER L(N)
C
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
