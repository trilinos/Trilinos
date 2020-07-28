C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK, FLGRAD
      SUBROUTINE FLGRAD(IEL,ICOUNT,IELLST,CNTRA,SHLNRM,SOLEA,SOLGRA,
     &                  ITT,iblk)

C  *********************************************************************

C  Subroutine FLGRAD performs the actual computation of the gradient
C  coefficients using the stuff that was set up in ELGRAD.
C  Start by translating into isoparametric coords (very helpful for
C  shells and makes some sense for continuum and is very useful for
C  subsequent interpolation.
C  Then do constrained least sqares to f=a0+a1eta+a2ksi+a3phi to
C  compute a1, a2, and a3 and stuff results into SOLGRA

C  Calls subroutines ERROR

C  Called by ELGRAD

C  *********************************************************************

C   IEL       the current element being worked
C   ICOUNT    the number of elements that share a node with IEL
C   IELLST    local list of elements that share a node with element
C             currently being processed
C   CNTRA     a list of element centroid coordinates for all elements
C             in the current element block (1:ndima,1:numeba)
C   SHLNRM    normal vector for shell element (1:3)
C   SOLEA     element variables (1:numeba,1:nvarel)
C   SOLGRA    element variable gradient (1:ndima,1:numeba,1:nvarel)
C   ITT       truth table
C   iblk      element block being processed (not ID)
C   IRED      reduction in dimensionality of data
C             IRED=0 full, IRED=1 colinear, IRED=2 coplanar

C  *********************************************************************

      include 'aexds1.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'

      DIMENSION CNTRA(NUMEBA,*), SHLNRM(3), SOLEA(NUMEBA,*)
      DIMENSION SOLGRA(NDIMA,NUMEBA,*), IELLST(100), ITT(NVAREL,*)
      DOUBLE PRECISION S(3,3), G(3), F(3), X(3)
      INTEGER L(3)

C  *********************************************************************

      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
        CALL ERROR('FLGRAD','ELEMENT TYPE',' ',ITYPE,
     &             'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &              0,' ',' ',1)
      END IF
      IRED = 0
      DO I = 1,3
         DO J = 1,3
            S(I,J) = 0.0D0
         end do
      end do

C  *********************************************************************
      IF (ITYPE .EQ. 13)THEN

C Shell element processing (quasi 2-D)

C If no elements connected, there can be no gradient

        IF (ICOUNT .EQ. 0)THEN
          DO 110 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 110
            SOLGRA(1,IEL,IVAR) = 0.
            SOLGRA(2,IEL,IVAR) = 0.
            SOLGRA(3,IEL,IVAR) = 0.
 110      CONTINUE
          GO TO 500

C If only one element connected, data is colinear
C If two elements connected, data may be colinear

        ELSE IF (ICOUNT .EQ. 1)THEN
          IRED = 1
        ELSE IF (ICOUNT .EQ. 2)THEN

C Check for colinearity. Create unit vector to 1st connected element
C centroid. Create unit vector to 2nd connected element centroid. Dot
C 1st unit vector with 2nd unit vector. If dot product is greater than
C 0.9, then data is colinear (IRED = 1)

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          V1MAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V13 = V13 / V1MAG
          V21 = CNTRA(IELLST(2),1) - CNTRA(IEL,1)
          V22 = CNTRA(IELLST(2),2) - CNTRA(IEL,2)
          V23 = CNTRA(IELLST(2),3) - CNTRA(IEL,3)
          V2MAG = SQRT(V21*V21 + V22*V22 + V23*V23)
          V21 = V21 / V2MAG
          V22 = V22 / V2MAG
          V23 = V23 / V2MAG
          VDOT = ABS(V11*V21 + V12*V22 + V13*V23)
          IF (VDOT .GT. 0.9)THEN
            IRED = 1
          END IF
        END IF
        IF (IRED .EQ. 1)THEN
C Colinear data

C rotate into vector
C NOTE: for colinearity, the XD coord is also the magnitude

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          V1MAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V13 = V13 / V1MAG

          S1 = 0.
          DO 120 I = 1, ICOUNT
            XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            XD = XORI*V11 + YORI*V12 +ZORI*V13
            S1 = S1 + (XD*XD)
 120      CONTINUE

          DO 130 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 130
            F1 = 0.
            DO 140 I = 1, ICOUNT
              XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
              XD = XORI*V11 + YORI*V12 +ZORI*V13
              F1 = SOLEA(IEL,IVAR) * XD
 140        CONTINUE
            X1D = F1 / S1
            SOLGRA(1,IEL,IVAR) = X1D * V11
            SOLGRA(2,IEL,IVAR) = X1D * V12
            SOLGRA(3,IEL,IVAR) = X1D * V13
 130      CONTINUE
        ELSE

C Rotate into plane of element and treat as if 2-D

C first unit vector

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          VMAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / VMAG
          V12 = V12 / VMAG
          V13 = V13 / VMAG

C compute 2nd (orthogonal) vector in plane - make it a unit vector

          V21 = V12 * SHLNRM(3) - SHLNRM(2) * V13
          V22 = SHLNRM(1) * V13 - V11 * SHLNRM(3)
          V23 = V11 * SHLNRM(2) - SHLNRM(1) * V12
          VMAG2 = SQRT(V21*V21 + V22*V22 + V23*V23)
          V21 = V21 / VMAG2
          V22 = V22 / VMAG2
          V23 = V23 / VMAG2

C Fill up matrix for constrained linear least squares

          DO 150 I = 1, ICOUNT

C rotate coords

            XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            XD = XORI*V11 + YORI*V12 + ZORI*V13
            YD = XORI*V21 + YORI*V22 + ZORI*V23

            S(1,1) = S(1,1) + DBLE(XD * XD)
            S(1,2) = S(1,2) + DBLE(XD * YD)
            S(2,2) = S(2,2) + DBLE(YD * YD)
 150      CONTINUE
          S(2,1) = S(1,2)
          S(1,3) = 0.D+00
          S(2,3) = 0.D+00
          S(3,3) = 1.D+00
          S(3,2) = 0.D+00
          S(3,1) = 0.D+00

C Forward Gauss

          CALL FRGE(3,S,L,G)

C Set up load vectors

          DO 160 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 160
            F(1) = 0.D+00
            F(2) = 0.D+00
            F(3) = 0.D+00
            DO 170 I = 1, ICOUNT
              XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
              XD = XORI*V11 + YORI*V12 + ZORI*V13
              YD = XORI*V21 + YORI*V22 + ZORI*V23

              F(1) = F(1) + DBLE(XD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
              F(2) = F(2) + DBLE(YD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
 170        CONTINUE

C Back substitution and load gradient coefficients into SOLGRA

            CALL BS(3,S,F,L,X)

C Rotate back and fill up gradient array

            SOLGRA(1,IEL,IVAR) = SNGL(X(2))*V11 + SNGL(X(3))*V21
            SOLGRA(2,IEL,IVAR) = SNGL(X(2))*V12 + SNGL(X(3))*V22
            SOLGRA(3,IEL,IVAR) = SNGL(X(2))*V13 + SNGL(X(3))*V23
 160      CONTINUE
        END IF

C  *********************************************************************
      ELSE IF (ITYPE .EQ. 3)THEN

C Quad element processing (2-D)

C If no elements connected, there can be no gradient

        IF (ICOUNT .EQ. 0)THEN
          DO 210 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 210
            SOLGRA(1,IEL,IVAR) = 0.
            SOLGRA(2,IEL,IVAR) = 0.
 210      CONTINUE
          GO TO 500

C If only one element connected, data is colinear
C If two elements connected, data may be colinear

        ELSE IF (ICOUNT .EQ. 1)THEN
          IRED = 1
        ELSE IF (ICOUNT .EQ. 2)THEN

C Check for colinearity. Create unit vector to 1st connected element
C centroid. Create unit vector to 2nd connected element centroid. Dot 1st
C unit vector with 2nd unit vector. If mag of dot product is greater than
C 0.9, then data is colinear (IRED = 1)

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V1MAG = SQRT(V11*V11 + V12*V12)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V21 = CNTRA(IELLST(2),1) - CNTRA(IEL,1)
          V22 = CNTRA(IELLST(2),2) - CNTRA(IEL,2)
          V2MAG = SQRT(V21*V21 + V22*V22)
          V21 = V21 / V2MAG
          V22 = V22 / V2MAG
          VDOT = ABS(V11*V21 + V12*V22)
          IF (VDOT .GT. 0.9)THEN
            IRED = 1
          END IF
        END IF
        IF(IRED .EQ. 1)THEN

C Colinear data
C Note: constraint at X=0 implies a0 so only 1 equation
C       remains to solve for a1

C rotate into vector
C NOTE: for colinearity, the XD coord is also the magnitude

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          VMAG = SQRT(V11*V11 + V12*V12)
          V11 = V11 / VMAG
          V12 = V12 / VMAG

          S1 = 0.
          DO 220 I = 1, ICOUNT
            XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            XD = XORI*V11 + YORI*V12
            S1 = S1 + (XD*XD)
 220      CONTINUE

          DO 230 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 230
            F1 = 0.
            DO 240 I = 1, ICOUNT
              XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              XD = XORI*V11 + YORI*V12
              F1 = SOLEA(IEL,IVAR) * XD
 240        CONTINUE
            X1D = F1 / S1
            SOLGRA(1,IEL,IVAR) = X1D * V11
            SOLGRA(2,IEL,IVAR) = X1D * V12
 230      CONTINUE
        ELSE IF(IRED .EQ. 0)THEN

C Full 2-D data

C Fill up matrix for constrained linear least squares

          DO 250 I = 1, ICOUNT
            XD = (CNTRA(IELLST(I),1) - CNTRA(IEL,1))
            YD = (CNTRA(IELLST(I),2) - CNTRA(IEL,2))
            S(1,1) = S(1,1) + DBLE(XD * XD)
            S(1,2) = S(1,2) + DBLE(XD * YD)
            S(2,2) = S(2,2) + DBLE(YD * YD)
 250      CONTINUE
          S(2,1) = S(1,2)
          S(1,3) = 0.D+00
          S(2,3) = 0.D+00
          S(3,3) = 1.D+00
          S(3,2) = 0.D+00
          S(3,1) = 0.D+00

C Forward Gauss

          CALL FRGE(3,S,L,G)

C Set up load vectors

          DO 260 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 260
            F(1) = 0.D+00
            F(2) = 0.D+00
            F(3) = 0.D+00
            DO 270 I = 1, ICOUNT
              XD = (CNTRA(IELLST(I),1) - CNTRA(IEL,1))
              YD = (CNTRA(IELLST(I),2) - CNTRA(IEL,2))
              F(1) = F(1) + DBLE(XD * (SOLEA(IELLST(I),IVAR)
     &                                - SOLEA(IEL,IVAR)))
              F(2) = F(2) + DBLE(YD * (SOLEA(IELLST(I),IVAR)
     &                                - SOLEA(IEL,IVAR)))
 270        CONTINUE

C Back substitution and load gradient coefficients into SOLGRA

            CALL BS(3,S,F,L,X)

C Fill up gradient array

            SOLGRA(1,IEL,IVAR) = SNGL(X(1))
            SOLGRA(2,IEL,IVAR) = SNGL(X(2))
 260      CONTINUE
        END IF

C  *********************************************************************
      ELSE IF (ITYPE .EQ. 10 .OR. ITYPE .EQ. 6)THEN

C Hex element processing (3-D)

C If no elements connected, there can be no gradient

        IF (ICOUNT .EQ. 0)THEN
          DO 310 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 310
            SOLGRA(1,IEL,IVAR) = 0.
            SOLGRA(2,IEL,IVAR) = 0.
            SOLGRA(3,IEL,IVAR) = 0.
 310      CONTINUE
          GO TO 500

C If only one element connected, data is colinear
C If two elements connected, data may be colinear

        ELSE IF (ICOUNT .EQ. 1)THEN
          IRED = 1
        ELSE IF (ICOUNT .EQ. 2)THEN

C Check for colinearity. Create unit vector to 1st connected element
C centroid. Create unit vector to 2nd connected element centroid. Dot 1st
C unit vector with 2nd unit vector. If mag of dot product is greater than
C 0.9, then data is colinear (IRED = 1)

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          V1MAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V13 = V13 / V1MAG
          V21 = CNTRA(IELLST(2),1) - CNTRA(IEL,1)
          V22 = CNTRA(IELLST(2),2) - CNTRA(IEL,2)
          V23 = CNTRA(IELLST(2),3) - CNTRA(IEL,3)
          V2MAG = SQRT(V21*V21 + V22*V22 + V23*V23)
          V21 = V21 / V2MAG
          V22 = V22 / V2MAG
          V23 = V23 / V2MAG
          VDOT = ABS(V11*V21 + V12*V22 + V13*V23)
          IF (VDOT .GT. 0.9)THEN
            IRED = 1
          ELSE
            IRED = 2

C  X-product vector-1 with vector-2 to get normal to plane defined
C  by the two vectors; then make it a unit vector

          VN1 = V12*V23 - V22*V13
          VN2 = V13*V21 - V11*V23
          VN3 = V11*V22 - V21*V12
          VNMAG = SQRT(VN1*VN1 + VN2*VN2 + VN3*VN3)
          END IF
        ELSE IF (ICOUNT .GT. 2 .AND. ICOUNT .LT. 11)THEN

C Check for coplanarity
C vector to 1st connected element

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)

C  Construct a vector to second connected element

          ICK = 2
 311      CONTINUE
          V21 = CNTRA(IELLST(ICK),1) - CNTRA(IEL,1)
          V22 = CNTRA(IELLST(ICK),2) - CNTRA(IEL,2)
          V23 = CNTRA(IELLST(ICK),3) - CNTRA(IEL,3)

C  X-product vector-1 with vector-2 to get normal to plane defined
C  by the two vectors; then make it a unit vector

          VN1 = V12*V23 - V22*V13
          VN2 = V13*V21 - V11*V23
          VN3 = V11*V22 - V21*V12
          VNMAG = SQRT(VN1*VN1 + VN2*VN2 + VN3*VN3)

C check for colinearity of elements
C if colinear, get new element and try again

          IF (VNMAG .LT. 1.E-13)THEN
            ICK = 3
            GO TO 311
          END IF
          VN1 = VN1 / VNMAG
          VN2 = VN2 / VNMAG
          VN3 = VN3 / VNMAG

C  Dot product of normal vector with unit vectors
C  to the remaining element centroids. If dot product
C  is too small, set IRED=2 and try the next vector.
C  If dot product is more than 0.1, data is not coplanar
C  set IRED = 0 and get on with it.

          DO 320 I = 3, ICOUNT
            V1 = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            V2 = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            V3 = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            VIMAG = SQRT(V1*V1 + V2*V2 + V3*V3)
            V1 = V1 / VIMAG
            V2 = V2 / VIMAG
            V3 = V3 / VIMAG
            VDOT = ABS(VN1*V1 + VN2*V2 + VN3*V3)
            IF (VDOT .LT. 0.1)THEN
              IRED = 2
            ELSE
              IRED = 0
              GO TO 330
            END IF
 320      CONTINUE

 330    CONTINUE
        END IF

        IF (IRED .EQ. 1)THEN

C Colinear data

C rotate into vector
C NOTE: for colinearity, the XD coord is also the magnitude

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          V1MAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V13 = V13 / V1MAG

          S1 = 0.
          DO 340 I = 1, ICOUNT
            XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            XD = XORI*V11 + YORI*V12 +ZORI*V13
            S1 = S1 + (XD*XD)
 340      CONTINUE

          DO 350 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 350
            F1 = 0.
            DO 360 I = 1, ICOUNT
              XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
              XD = XORI*V11 + YORI*V12 +ZORI*V13
              F1 = SOLEA(IEL,IVAR) * XD
 360        CONTINUE
            X1D = F1 / S1
            SOLGRA(1,IEL,IVAR) = X1D * V11
            SOLGRA(2,IEL,IVAR) = X1D * V12
            SOLGRA(3,IEL,IVAR) = X1D * V13
 350      CONTINUE

        ELSE IF (IRED .EQ. 2)THEN

C Coplanar data

C first unit vector

          V11 = CNTRA(IELLST(1),1) - CNTRA(IEL,1)
          V12 = CNTRA(IELLST(1),2) - CNTRA(IEL,2)
          V13 = CNTRA(IELLST(1),3) - CNTRA(IEL,3)
          V1MAG = SQRT(V11*V11 + V12*V12 + V13*V13)
          V11 = V11 / V1MAG
          V12 = V12 / V1MAG
          V13 = V13 / V1MAG

C compute 2nd (orthogonal) vector in plane - make it a unit vector

          V21 = V12 * VN3 - VN2 * V13
          V22 = VN1 * V13 - V11 * VN3
          V23 = V11 * VN2 - VN1 * V12
          V2MAG = SQRT(V21*V21 + V22*V22 + V23*V23)
          V21 = V21 / V2MAG
          V22 = V22 / V2MAG
          V23 = V23 / V2MAG

C Fill up matrix for constrained linear least squares

          DO 370 I = 1, ICOUNT

C rotate coords

            XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            XD = XORI*V11 + YORI*V12 + ZORI*V13
            YD = XORI*V21 + YORI*V22 + ZORI*V23

            S(1,1) = S(1,1) + DBLE(XD * XD)
            S(1,2) = S(1,2) + DBLE(XD * YD)
            S(2,2) = S(2,2) + DBLE(YD * YD)
 370      CONTINUE
          S(2,1) = S(1,2)
          S(1,3) = 0.D+00
          S(2,3) = 0.D+00
          S(3,3) = 1.D+00
          S(3,2) = 0.D+00
          S(3,1) = 0.D+00

C Forward Gauss

          CALL FRGE(3,S,L,G)

C Set up load vectors

          DO 380 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 380
            F(1) = 0.D+00
            F(2) = 0.D+00
            F(3) = 0.D+00
            DO 390 I = 1, ICOUNT
              XORI = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YORI = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              ZORI = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
              XD = XORI*V11 + YORI*V12 + ZORI*V13
              YD = XORI*V21 + YORI*V22 + ZORI*V23

              F(1) = F(1) + DBLE(XD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
              F(2) = F(2) + DBLE(YD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
 390        CONTINUE

C Back substitution and load gradient coefficients into SOLGRA

            CALL BS(3,S,F,L,X)

C Rotate back and fill up gradient array

            SOLGRA(1,IEL,IVAR) = SNGL(X(1))*V11 + SNGL(X(2))*V21
            SOLGRA(2,IEL,IVAR) = SNGL(X(1))*V12 + SNGL(X(2))*V22
            SOLGRA(3,IEL,IVAR) = SNGL(X(1))*V13 + SNGL(X(2))*V23
 380      CONTINUE

        ELSE

C Fully 3-D data (IRED=0)

C Zero the matrix
C Note: constraint at X=Y=Z=0 implies a0 so only 3 equations
C       remain to solve for a1, a2, and a3

          DO I = 1, 3
             DO J = 1, 3
                S(I,J) = 0.D+00
             end do
          end do

C Fill up matrix for constrained linear least squares

          DO 410 I = 1, ICOUNT
            XD = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
            YD = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
            ZD = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
            S(1,1) = S(1,1) + DBLE(XD * XD)
            S(1,2) = S(1,2) + DBLE(XD * YD)
            S(1,3) = S(1,3) + DBLE(XD * ZD)
            S(2,2) = S(2,2) + DBLE(YD * YD)
            S(2,3) = S(2,3) + DBLE(YD * ZD)
            S(3,3) = S(3,3) + DBLE(ZD * ZD)
 410      CONTINUE
          S(2,1) = S(1,2)
          S(3,1) = S(1,3)
          S(3,2) = S(2,3)

C Forward Gauss

          CALL FRGE(3,S,L,G)

C Set up load vectors

          DO 420 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 420
            F(1) = 0.D+00
            F(2) = 0.D+00
            F(3) = 0.D+00
            DO 430 I = 1, ICOUNT
              XD = CNTRA(IELLST(I),1) - CNTRA(IEL,1)
              YD = CNTRA(IELLST(I),2) - CNTRA(IEL,2)
              ZD = CNTRA(IELLST(I),3) - CNTRA(IEL,3)
              F(1) = F(1) + DBLE(XD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
              F(2) = F(2) + DBLE(YD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
              F(3) = F(3) + DBLE(ZD * (SOLEA(IELLST(I),IVAR)
     &                                 - SOLEA(IEL,IVAR)))
 430        CONTINUE

C Back substitution and load gradient coefficients into SOLGRA

            CALL BS(3,S,F,L,X)

C Fill up gradient array

            SOLGRA(1,IEL,IVAR) = SNGL(X(1))
            SOLGRA(2,IEL,IVAR) = SNGL(X(2))
            SOLGRA(3,IEL,IVAR) = SNGL(X(3))
 420      CONTINUE
        END IF
      ELSE
        CALL ERROR ('MAPVAR','INCORRECT ELEMENT TYPE',
     &              'ELEMENT TYPE =',ITYPE,
     &              'NOT YET IMPLEMENTED',0,' ',' ',1)
      END IF
 500  CONTINUE
      RETURN
      END
