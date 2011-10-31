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
      SUBROUTINE EXTH(IGLND,INVCN,MAXLN,NOD,INVLEN,XA,YA,ZA,CNTRA,
     &                SOLEA,SOLENA,ITT,iblk)
C
C************************************************************************
C
C Subroutine EXTH sets up the matrix and vectors for a least squares
C linear interpolation/extrapolation of element variable data to the 
C nodes for 3-D elements. This routine has been checked out for 8-node
C hex and 4-node and 8-node (treated same as 4-node) tet elements.
C In the special case of data from only 4 elements, the result is not 
C a true least squares fit in that the least squares error is zero.
C
C
C Calls subroutines FRGE & BS
C
C Called by ELTON3
C
C************************************************************************
C
C  IGLND  INT   The global node number being processed
C  INVCN  INT   Inverse connectivity (1:maxln,1:numnda)
C  MAXLN  INT   The maximum number of elements connected to any node
C  NOD    INT   The local node used to get elements from INVCN
C  INVLEN INT   The number of elements connected to NOD
C  XA,etc REAL  Vectors containing nodal coordinates
C  CNTRA  REAL  Array containing the coordinates of the element 
C               centroids (1:3)
C  SOLEA  REAL  The element variables
C  SOLENA REAL  Element variables at nodes
C               number with the global mesh node number (1:numnda)
C  ITT    INT   truth table
C  iblk   INT   element block being processed (not ID)
C  ICOP   INT   Flag defining co-planarity of elements (0-coplan, 1-not)
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
      DIMENSION SOLENA(NODESA,NVAREL),ITT(NVAREL,*)
      DOUBLE PRECISION S(4,4),G(4),F(4),X(4)
      INTEGER L(4)
C
C************************************************************************
      ICOP = 0
C
C  First check elements for coplanarity
C
C  Construct a vector from first element centroid to second
C
      VEC11 = CNTRA(INVCN(2,NOD),1) - CNTRA(INVCN(1,NOD),1)
      VEC12 = CNTRA(INVCN(2,NOD),2) - CNTRA(INVCN(1,NOD),2)
      VEC13 = CNTRA(INVCN(2,NOD),3) - CNTRA(INVCN(1,NOD),3)
      V1MAG = SQRT(VEC11*VEC11 + VEC12*VEC12 + VEC13*VEC13)
C
C  Construct a vector from first element centroid to third
C
      VEC21 = CNTRA(INVCN(3,NOD),1) - CNTRA(INVCN(1,NOD),1)
      VEC22 = CNTRA(INVCN(3,NOD),2) - CNTRA(INVCN(1,NOD),2)
      VEC23 = CNTRA(INVCN(3,NOD),3) - CNTRA(INVCN(1,NOD),3)
C
C  X-product vector-1 with vector-2 to get normal to plane defined
C  by the two vectors then make a unit vector
C
      VN1 = VEC12*VEC23 - VEC22*VEC13
      VN2 = VEC13*VEC21 - VEC11*VEC23
      VN3 = VEC11*VEC22 - VEC21*VEC12
      VNMAG = SQRT(VN1*VN1 + VN2*VN2 + VN3*VN3)
      VN1 = VN1 / VNMAG
      VN2 = VN2 / VNMAG
      VN3 = VN3 / VNMAG
C
C  Dot product of normal vector with vectors from first element 
C  centroid to the remaining element centroids. If dot product
C  is too small, data is coplanar - try the next vector.
C  If dot product more then 0.1 times the vector, data is not
C  coplanar, set ICOP to 1 and get on with it.
C
      DO 5 I = 4, INVLEN
        VEC1 = CNTRA(INVCN(I,NOD),1) - CNTRA(INVCN(1,NOD),1)
        VEC2 = CNTRA(INVCN(I,NOD),2) - CNTRA(INVCN(1,NOD),2)
        VEC3 = CNTRA(INVCN(I,NOD),3) - CNTRA(INVCN(1,NOD),3)
        VIMAG = SQRT(VEC1*VEC1 + VEC2*VEC2 + VEC3*VEC3)
        VDTMAG = SQRT(VN1*VEC1*VN1*VEC1 + VN2*VEC2*VN2*VEC2
     &              + VN3*VEC3*VN3*VEC3)
        COMP = VDTMAG * 10.
        IF (COMP .LT. VIMAG)THEN
          GO TO 5
        ELSE
          ICOP = 1
          go to 6
        END IF
 5    CONTINUE
C
C  Zero matrix
C
 6    CONTINUE
      DO 10 I = 1,4
      DO 10 J = 1,4
        S(I,J) = 0.D+00
   10 CONTINUE
C
C Branch on coplanar data vs truely 3-d data
C
      IF (ICOP .EQ. 1)THEN
C
C  Set up matrix for linear fit
C
        S(1,1) = DBLE(INVLEN)
        DO 20 I = 1, INVLEN
          S(1,2) = S(1,2)+DBLE(XA(IGLND) - CNTRA(INVCN(I,NOD),1))
          S(1,3) = S(1,3)+DBLE(YA(IGLND) - CNTRA(INVCN(I,NOD),2))
          S(1,4) = S(1,4)+DBLE(ZA(IGLND) - CNTRA(INVCN(I,NOD),3))
          S(2,2) = S(2,2)+DBLE((XA(IGLND) - CNTRA(INVCN(I,NOD),1)) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
          S(2,3) = S(2,3)+DBLE((YA(IGLND) - CNTRA(INVCN(I,NOD),2)) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
          S(2,4) = S(2,4)+DBLE((ZA(IGLND) - CNTRA(INVCN(I,NOD),3)) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
          S(3,3) = S(3,3)+DBLE((YA(IGLND) - CNTRA(INVCN(I,NOD),2)) *
     &                    (YA(IGLND) - CNTRA(INVCN(I,NOD),2)))
          S(3,4) = S(3,4)+DBLE((ZA(IGLND) - CNTRA(INVCN(I,NOD),3)) *
     &                    (YA(IGLND) - CNTRA(INVCN(I,NOD),2)))
          S(4,4) = S(4,4)+DBLE((ZA(IGLND) - CNTRA(INVCN(I,NOD),3)) *
     &                    (ZA(IGLND) - CNTRA(INVCN(I,NOD),3)))
   20   CONTINUE
        S(2,1) = S(1,2)
        S(3,1) = S(1,3)
        S(4,1) = S(1,4)
        S(3,2) = S(2,3)
        S(4,2) = S(2,4)
        S(4,3) = S(3,4)
C
C  Forward Gauss elimination (Kincaid pg. 220) (double precision)
C
        CALL FRGE(4,S,L,G)
C
C  Set up load vectors - number of element variables
C
        DO 30 IVAR = 1, NVAREL
          IF (ITT(IVAR,iblk) .EQ. 0)GO TO 30
          F(1) = 0.D+00
          F(2) = 0.D+00
          F(3) = 0.D+00
          F(4) = 0.D+00
          DO 40 I = 1, INVLEN
            F(1) = F(1) + DBLE(SOLEA(INVCN(I,NOD),IVAR))
            F(2) = F(2) + DBLE(SOLEA(INVCN(I,NOD),IVAR) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
            F(3) = F(3) + DBLE(SOLEA(INVCN(I,NOD),IVAR) *
     &                    (YA(IGLND) - CNTRA(INVCN(I,NOD),2)))
            F(4) = F(4) + DBLE(SOLEA(INVCN(I,NOD),IVAR) *
     &                    (ZA(IGLND) - CNTRA(INVCN(I,NOD),3)))
   40     CONTINUE
C
C  Back substitution (Kincaid pg. 223) (double precision)
C
          CALL BS(4,S,F,L,X)
C
C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated to (IGLND), thus X, Y, Z are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y + X(4) * Z
C
          SOLENA(IGLND,IVAR) = SNGL(X(1))
   30   CONTINUE
C
      ELSE IF (ICOP .EQ. 0)THEN
C
C first unit vector
C
        V11 = VEC11 / V1MAG
        V12 = VEC12 / V1MAG
        V13 = VEC13 / V1MAG
C
C compute 2nd (orthogonal) vector in plane - make it a unit vector
C
        V21 = V12 * VN3 - VN2 * V13
        V22 = VN1 * V13 - V11 * VN3
        V23 = V11 * VN2 - VN1 * V12
        V2MAG = SQRT(V21*V21 + V22*V22 + V23*V23)
        V21 = V21 / V2MAG
        V22 = V22 / V2MAG
        V23 = V23 / V2MAG
C
C set up matrix for least squares fit
C
        S(1,1) = DBLE(INVLEN)
        DO 50 I = 1, INVLEN
C
C rotate coords
C
          XORI = XA(IGLND)-CNTRA(INVCN(I,NOD),1)
          YORI = YA(IGLND)-CNTRA(INVCN(I,NOD),2)
          ZORI = ZA(IGLND)-CNTRA(INVCN(I,NOD),3)
          XP = XORI*V11 + YORI*V12 + ZORI*V13
          YP = XORI*V21 + YORI*V22 + ZORI*V23
C
          S(1,2) = S(1,2)+DBLE(XP)
          S(1,3) = S(1,3)+DBLE(YP)
          S(2,2) = S(2,2)+DBLE(XP * XP)
          S(2,3) = S(2,3)+DBLE(XP * YP)
          S(3,3) = S(3,3)+DBLE(YP * YP)
 50     CONTINUE
        S(2,1) = S(1,2)
        S(3,1) = S(1,3)
        S(3,2) = S(2,3)
        S(1,4) = 0.D+00
        S(2,4) = 0.D+00
        S(3,4) = 0.D+00
        S(4,4) = 1.D+00
        S(4,3) = 0.D+00
        S(4,2) = 0.D+00
        S(4,1) = 0.D+00
C
C  Forward Gauss elimination (Kincaid pg. 220) (double precision)
C
        CALL FRGE(4,S,L,G)
C
C  Set up load vectors - number of element variables
C
        DO 60 IVAR = 1, NVAREL
          IF (ITT(IVAR,iblk) .EQ. 0)GO TO 60
          F(1) = 0.D+00
          F(2) = 0.D+00
          F(3) = 0.D+00
          F(4) = 0.D+00
          DO 70 I = 1, INVLEN
            XORI = XA(IGLND)-CNTRA(INVCN(I,NOD),1)
            YORI = YA(IGLND)-CNTRA(INVCN(I,NOD),2)
            ZORI = ZA(IGLND)-CNTRA(INVCN(I,NOD),3)
            XP = XORI*V11 + YORI*V12 + ZORI*V13
            YP = XORI*V21 + YORI*V22 + ZORI*V23
            F(1) = F(1) + DBLE(SOLEA(INVCN(I,NOD),IVAR))
            F(2) = F(2) + DBLE(SOLEA(INVCN(I,NOD),IVAR) * XP)
            F(3) = F(3) + DBLE(SOLEA(INVCN(I,NOD),IVAR) * YP)
   70     CONTINUE
C
C  Back substitution (Kincaid pg. 223) (double precision)
C
          CALL BS(4,S,F,L,X)
C
C Ordinaly you would need to rotate back into cartesian coords
C however, we only need X(1) so there is no need to rotate here
C
C          X2 = SNGL(X(2))*V11 + SNGL(X(3))*V21
C          X3 = SNGL(X(2))*V12 + SNGL(X(3))*V22
C          X4 = SNGL(X(2))*V13 + SNGL(X(3))*V23
C          X(2) = X2
C          X(3) = X3
C          X(4) = X4
C
C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated to (IGLND), thus X, Y, Z are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y + X(4) * Z
C
          SOLENA(IGLND,IVAR) = SNGL(X(1))
   60   CONTINUE
      END IF
      RETURN
      END
