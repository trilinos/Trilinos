C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C========================================================================
      SUBROUTINE EXTQ(IGLND,INVCN,MAXLN,NOD,INVLEN,XA,YA,CNTRA,SOLEA,
     &                SOLENA,ITT,iblk)

C************************************************************************

C Subroutine EXTQ sets up the matrix and vectors for a least squares
C linear interpolation/extrapolation of element variable data to the
C nodes for a 4-node quad element. In the special case of data from
C only 3 elements, the result is not least squares fit but a
C triangularization.

C Calls subroutines FRGE & BS

C Called by ELTON3

C************************************************************************

C  IGLND  INT   The global node number being processed
C  INVCN  INT   The inverse connectivity (1:maxln,1:numnda)
C  MAXLN  INT   The maximum number of elements connected to any node
C  NOD    INT   The node used to get the elements from INVCN
C  INVLEN INT   The number of elements connected to NOD
C  XA,etc REAL  Vectors containing nodal coordinates
C  CNTRA  REAL  Array containing the coordinates of the element
C               centroids (1:3)
C  SOLEA  REAL  The element variables
C  SOLENA REAL  Element variables at nodes
C               number with the global mesh node number (1:numnda)
C  ITT    INT   truth table
C  IM     INT   element block being processed (not ID)
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

      DIMENSION INVCN(MAXLN,*),XA(*),YA(*)
      DIMENSION CNTRA(NUMEBA,*),SOLEA(NUMEBA,*)
      DIMENSION SOLENA(NODESA,NVAREL), ITT(NVAREL,*)
      DOUBLE PRECISION S(3,3),G(3),F(3),X(3)
      INTEGER L(3)

C************************************************************************

C  Zero matrix

      DO I = 1,3
         DO J = 1,3
            S(I,J) = 0.D+00
         end do
      end do

C  Set up matrix for linear fit

      S(1,1) = DBLE(INVLEN)
      DO I = 1, INVLEN
        S(1,2) = S(1,2) + DBLE(XA(IGLND) - CNTRA(INVCN(I,NOD),1))
        S(1,3) = S(1,3) + DBLE(YA(IGLND) - CNTRA(INVCN(I,NOD),2))
        S(2,2) = S(2,2) + DBLE((XA(IGLND) - CNTRA(INVCN(I,NOD),1)) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
        S(2,3) = S(2,3) + DBLE((YA(IGLND) - CNTRA(INVCN(I,NOD),2)) *
     &                    (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
        S(3,3) = S(3,3) + DBLE((YA(IGLND) - CNTRA(INVCN(I,NOD),2)) *
     &                    (YA(IGLND) - CNTRA(INVCN(I,NOD),2)))
      end do
      S(2,1) = S(1,2)
      S(3,1) = S(1,3)
      S(3,2) = S(2,3)

C  Forward Gauss elimination (Kincaid pg. 220) (double precision)

      CALL FRGE(3,S,L,G)

C  Set up load vectors - number of element variables

      DO IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 30
        F(1) = 0.D+00
        F(2) = 0.D+00
        F(3) = 0.D+00
        DO I = 1, INVLEN
          F(1) = F(1) + DBLE(SOLEA(INVCN(I,NOD),IVAR))
          F(2) = F(2) + DBLE(SOLEA(INVCN(I,NOD),IVAR) *
     &                  (XA(IGLND) - CNTRA(INVCN(I,NOD),1)))
          F(3) = F(3) + DBLE(SOLEA(INVCN(I,NOD),IVAR) *
     &                  (YA(IGLND) - CNTRA(INVCN(I,NOD),2)))
       end do

C  Back substitution (Kincaid pg. 223) (double precision)

        CALL BS(3,S,F,L,X)

C  Fill in nodal element value array (SOLENA)
C  Note: X and Y distances in S and F are centered on node being
C        interpolated to (IGLND), thus X and Y are zero in the eq.
C        Value = X(1) + X(2) * X + X(3) * Y

        SOLENA(IGLND,IVAR) = SNGL(X(1))
      end do
   30 CONTINUE
      RETURN
      END
