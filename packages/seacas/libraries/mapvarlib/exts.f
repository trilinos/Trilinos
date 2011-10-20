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

C========================================================================
      SUBROUTINE EXTS(IGLND,INVCN,MAXLN,NXGLND,INVLEN,XA,YA,ZA,
     &                CNTRA,SOLEA,SOLENA,ITT,iblk)
C
C************************************************************************
C
C Subroutine EXTS sets up the matrix and vectors for a least squares
C linear interpolation/extrapolation of element variable data to the 
C nodes for a 4-node quad element. In the special case of data from 
C only 3 elements, the result is not least squares fit but a 
C triangularization.
C
C Calls subroutines FRGE & BS
C
C Called by SELTN3
C
C************************************************************************
C
C  IGLND  INT   The global node number being processed
C  INVCN  INT   Inverse connectivity (1:maxln,1:numnda)
C  MAXLN  INT   The maximum number of elements connected to any node
C  NXGLND    INT   The local node used to get elements from INVCN
C  INVLEN INT   The number of elements connected to NXGLND
C  XA,etc REAL  Vectors containing nodal coordinates
C  CNTRA  REAL  Array containing the coordinates of the element 
C               centroids (1:3)
C  SOLEA  REAL  The element variables
C  SOLENA REAL  Element variables at nodes
C               number with the global mesh node number (1:numnda)
C  ITT    INT   truth table
C  iblk   INT   element block being processed (not ID)
C  INTND  INT   The global node number associated with IGLND
C  S      REAL  The coefficient matrix for the least squares fit
C  L      INT   Dummy vector - used in FRGE and BS
C  X      REAL  The solution vector - used in BS
C  G      REAL  Dummy vector - used in FRGE
C  F      REAL  The load vector for the least squares fit
C
C************************************************************************
C
      include 'aexds1.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'tapes.blk'
C
      DIMENSION INVCN(MAXLN,*),XA(*),YA(*),ZA(*)
      DIMENSION CNTRA(NUMEBA,*),SOLEA(NUMEBA,*)
      DIMENSION SOLENA(NODESA,NVAREL), ITT(NVAREL,*)
      DIMENSION IFRST(3), RLENTH(8), XLC(8), YLC(8)
C      DIMENSION ZLC(8)
      DOUBLE PRECISION S(3,3),G(3),F(3),X(3)
      INTEGER L(3)
C
C************************************************************************
C
C  Zero matrix
C
      DO 10 I = 1,3
        IFRST(I) = I
      DO 10 J = 1,3
        S(I,J) = 0.D+00
   10 CONTINUE
c
c find distance from interpolation point to element centroids
c
      DO 20 I = 1, INVLEN
        A = XA(IGLND) - CNTRA(INVCN(I,NXGLND),1)
        B = YA(IGLND) - CNTRA(INVCN(I,NXGLND),2)
        C = ZA(IGLND) - CNTRA(INVCN(I,NXGLND),3)
        RLENTH(I) = SQRT(A*A + B*B + C*C)
   20 CONTINUE
C
C find the three closest element centroids
C
      IF (INVLEN .EQ. 3) THEN
        DO 30 I = 1, 2
          IF (RLENTH(I) .GT. RLENTH(I+1))THEN
            ITEMP = IFRST(I)
            IFRST(I) = IFRST(I+1)
            IFRST(I+1) = ITEMP
          END IF
   30   CONTINUE
          IF (RLENTH(1) .GT. RLENTH(2))THEN
            ITEMP = IFRST(1)
            IFRST(1) = IFRST(2)
            IFRST(2) = ITEMP
          END IF
C
      ELSE
        DO 40 I = 2, INVLEN
          IF (RLENTH(I) .LT. RLENTH(IFRST(1))) IFRST(1) = I
   40   CONTINUE
        IFRST(2) = 1
        IF (IFRST(1) .EQ. 1) IFRST(2) = 2
        DO 50 I = IFRST(2), INVLEN
          IF (I .EQ. IFRST(1))GO TO 50
          IF (RLENTH(I) .LT. RLENTH(IFRST(2)))IFRST(2) = I
   50   CONTINUE
        IFRST(3) = 1
        IF (IFRST(1) .EQ. 1 .OR. IFRST(2) .EQ. 1) IFRST(3)=2
        IF (IFRST(1) .EQ. IFRST(3) .OR. IFRST(2) .EQ. IFRST(3))
     &  IFRST(3)=3
        DO 60 I = IFRST(3), INVLEN
          IF (I .EQ. IFRST(1))GO TO 60
          IF (I .EQ. IFRST(2))GO TO 60
          IF (RLENTH(I) .LT. RLENTH(IFRST(3))) IFRST(3) = I
   60   CONTINUE
      END IF
C
C use three closest element centroids to define a plane
C establish coordinate system on this plane centered on 
C interpolation point
C
      A11 = CNTRA(INVCN(IFRST(2),NXGLND),1) - 
     &      CNTRA(INVCN(IFRST(1),NXGLND),1)
      A12 = CNTRA(INVCN(IFRST(2),NXGLND),2) - 
     &      CNTRA(INVCN(IFRST(1),NXGLND),2)
      A13 = CNTRA(INVCN(IFRST(2),NXGLND),3) - 
     &      CNTRA(INVCN(IFRST(1),NXGLND),3)
      RLN = SQRT(A11*A11 + A12*A12 + A13*A13)
      A11 = A11/RLN
      A12 = A12/RLN
      A13 = A13/RLN
C
      A31 = (CNTRA(INVCN(IFRST(2),NXGLND),2) -
     &       CNTRA(INVCN(IFRST(1),NXGLND),2))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),3) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),3))
     &    - (CNTRA(INVCN(IFRST(2),NXGLND),3) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),3))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),2) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),2))
      A32 = (CNTRA(INVCN(IFRST(2),NXGLND),3) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),3))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),1) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),1))
     &    - (CNTRA(INVCN(IFRST(2),NXGLND),1) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),1))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),3) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),3))
      A33 = (CNTRA(INVCN(IFRST(2),NXGLND),1) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),1))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),2) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),2))
     &    - (CNTRA(INVCN(IFRST(2),NXGLND),2) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),2))
     &    * (CNTRA(INVCN(IFRST(3),NXGLND),1) - 
     &       CNTRA(INVCN(IFRST(1),NXGLND),1))
      RLN = SQRT(A31*A31 + A32*A32 + A33*A33)
      A31 = A31/RLN
      A32 = A32/RLN
      A33 = A33/RLN
C
      A21 = A32*A13 - A33*A12
      A22 = A11*A33 - A31*A13 
      A23 = A31*A12 - A11*A32
C
      DO 70 I = 1, INVLEN
          XLC(I) = A11 * (CNTRA(INVCN(I,NXGLND),1) - XA(IGLND))
     &           + A12 * (CNTRA(INVCN(I,NXGLND),2) - YA(IGLND))
     &           + A13 * (CNTRA(INVCN(I,NXGLND),3) - ZA(IGLND))
          YLC(I) = A21 * (CNTRA(INVCN(I,NXGLND),1) - XA(IGLND))
     &           + A22 * (CNTRA(INVCN(I,NXGLND),2) - YA(IGLND))
     &           + A23 * (CNTRA(INVCN(I,NXGLND),3) - ZA(IGLND))
c          ZLC(I) = A31 * (CNTRA(INVCN(I,NXGLND),1) - XA(IGLND))
c    &            + A32 * (CNTRA(INVCN(I,NXGLND),2) - YA(IGLND))
c    &            + A33 * (CNTRA(INVCN(I,NXGLND),3) - ZA(IGLND))
   70 CONTINUE
C
C
C  Set up matrix for linear fit
C
      S(1,1) = INVLEN
      DO 80 I = 1, INVLEN
        S(1,2) = S(1,2) + DBLE(XLC(I))
        S(1,3) = S(1,3) + DBLE(YLC(I))
        S(2,2) = S(2,2) + DBLE(XLC(I) * XLC(I))
        S(2,3) = S(2,3) + DBLE(YLC(I) * XLC(I))
        S(3,3) = S(3,3) + DBLE(YLC(I) * YLC(I))
   80 CONTINUE
      S(2,1) = S(1,2)
      S(3,1) = S(1,3)
      S(3,2) = S(2,3)
C
C  Forward Gauss elimination (Kincaid pg. 220) (double precision)
C
      CALL FRGE(3,S,L,G)
C
C  Set up load vectors - number of element variables
C
      DO 90 IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 90
        F(1) = 0.D+00
        F(2) = 0.D+00
        F(3) = 0.D+00
        DO 100 I = 1, INVLEN
          F(1) = F(1) + DBLE(SOLEA(INVCN(I,NXGLND),IVAR))
          F(2) = F(2) + DBLE(SOLEA(INVCN(I,NXGLND),IVAR) * XLC(I))
          F(3) = F(3) + DBLE(SOLEA(INVCN(I,NXGLND),IVAR) * YLC(I))
  100   CONTINUE
C
C  Back substitution (Kincaid pg. 223) (double precision)
C
        CALL BS(3,S,F,L,X)
C
C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated, thus X and Y are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y
C
        SOLENA(IGLND,IVAR) = SNGL(X(1))
   90 CONTINUE
      RETURN
      END
