C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C========================================================================
      SUBROUTINE EXTH(IGLND,INVCN,MAXLN,NOD,INVLEN,XA,YA,ZA,CNTRA,
     &                SOLEA,SOLENA,ITT,iblk)

C************************************************************************

C Subroutine EXTH sets up the matrix and vectors for a least squares
C linear interpolation/extrapolation of element variable data to the
C nodes for 3-D elements. This routine has been checked out for 8-node
C hex and 4-node and 8-node (treated same as 4-node) tet elements.
C In the special case of data from only 4 elements, the result is not
C a true least squares fit in that the least squares error is zero.

C Calls subroutines FRGE & BS

C Called by ELTON3

C************************************************************************

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

C************************************************************************

      include 'aexds1.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'tapes.blk'

      DIMENSION INVCN(MAXLN,*),XA(*),YA(*),ZA(*)
      DIMENSION CNTRA(NUMEBA,*),SOLEA(NUMEBA,*)
      DIMENSION SOLENA(NODESA,NVAREL),ITT(NVAREL,*)
      DOUBLE PRECISION S(4,4),G(4),F(4),X(4)
      INTEGER L(4)

C************************************************************************
      ICOP = 0

C  First check elements for coplanarity

C  Construct a vector from first element centroid to second

      VEC11 = CNTRA(INVCN(2,NOD),1) - CNTRA(INVCN(1,NOD),1)
      VEC12 = CNTRA(INVCN(2,NOD),2) - CNTRA(INVCN(1,NOD),2)
      VEC13 = CNTRA(INVCN(2,NOD),3) - CNTRA(INVCN(1,NOD),3)
      V1MAG = SQRT(VEC11*VEC11 + VEC12*VEC12 + VEC13*VEC13)

C  Construct a vector from first element centroid to third

      VEC21 = CNTRA(INVCN(3,NOD),1) - CNTRA(INVCN(1,NOD),1)
      VEC22 = CNTRA(INVCN(3,NOD),2) - CNTRA(INVCN(1,NOD),2)
      VEC23 = CNTRA(INVCN(3,NOD),3) - CNTRA(INVCN(1,NOD),3)

C  X-product vector-1 with vector-2 to get normal to plane defined
C  by the two vectors then make a unit vector

      VN1 = VEC12*VEC23 - VEC22*VEC13
      VN2 = VEC13*VEC21 - VEC11*VEC23
      VN3 = VEC11*VEC22 - VEC21*VEC12
      VNMAG = SQRT(VN1*VN1 + VN2*VN2 + VN3*VN3)
      VN1 = VN1 / VNMAG
      VN2 = VN2 / VNMAG
      VN3 = VN3 / VNMAG

C  Dot product of normal vector with vectors from first element
C  centroid to the remaining element centroids. If dot product
C  is too small, data is coplanar - try the next vector.
C  If dot product more then 0.1 times the vector, data is not
C  coplanar, set ICOP to 1 and get on with it.

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

C  Zero matrix

 6    CONTINUE
      DO I = 1,4
         DO J = 1,4
            S(I,J) = 0.D+00
         end do
      end do

C Branch on coplanar data vs truly 3-d data

      IF (ICOP .EQ. 1)THEN

C  Set up matrix for linear fit

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

C  Forward Gauss elimination (Kincaid pg. 220) (double precision)

        CALL FRGE(4,S,L,G)

C  Set up load vectors - number of element variables

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

C  Back substitution (Kincaid pg. 223) (double precision)

          CALL BS(4,S,F,L,X)

C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated to (IGLND), thus X, Y, Z are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y + X(4) * Z

          SOLENA(IGLND,IVAR) = SNGL(X(1))
   30   CONTINUE

      ELSE IF (ICOP .EQ. 0)THEN

C first unit vector

        V11 = VEC11 / V1MAG
        V12 = VEC12 / V1MAG
        V13 = VEC13 / V1MAG

C compute 2nd (orthogonal) vector in plane - make it a unit vector

        V21 = V12 * VN3 - VN2 * V13
        V22 = VN1 * V13 - V11 * VN3
        V23 = V11 * VN2 - VN1 * V12
        V2MAG = SQRT(V21*V21 + V22*V22 + V23*V23)
        V21 = V21 / V2MAG
        V22 = V22 / V2MAG
        V23 = V23 / V2MAG

C set up matrix for least squares fit

        S(1,1) = DBLE(INVLEN)
        DO 50 I = 1, INVLEN

C rotate coords

          XORI = XA(IGLND)-CNTRA(INVCN(I,NOD),1)
          YORI = YA(IGLND)-CNTRA(INVCN(I,NOD),2)
          ZORI = ZA(IGLND)-CNTRA(INVCN(I,NOD),3)
          XP = XORI*V11 + YORI*V12 + ZORI*V13
          YP = XORI*V21 + YORI*V22 + ZORI*V23

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

C  Forward Gauss elimination (Kincaid pg. 220) (double precision)

        CALL FRGE(4,S,L,G)

C  Set up load vectors - number of element variables

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

C  Back substitution (Kincaid pg. 223) (double precision)

          CALL BS(4,S,F,L,X)

C Ordinaly you would need to rotate back into cartesian coords
C however, we only need X(1) so there is no need to rotate here

C          X2 = SNGL(X(2))*V11 + SNGL(X(3))*V21
C          X3 = SNGL(X(2))*V12 + SNGL(X(3))*V22
C          X4 = SNGL(X(2))*V13 + SNGL(X(3))*V23
C          X(2) = X2
C          X(3) = X3
C          X(4) = X4

C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated to (IGLND), thus X, Y, Z are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y + X(4) * Z

          SOLENA(IGLND,IVAR) = SNGL(X(1))
   60   CONTINUE
      END IF
      RETURN
      END
