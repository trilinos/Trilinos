C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CUTROT (CUTPLA, CUTMAT, *)
C=======================================================================

C   --*** CUTROT *** (MESH) Get rotation matrix for cut
C   --   Written by Amy Gilkey - revised 07/03/86
C   --
C   --CUTROT makes up the rotation matrix that will rotate the mesh
C   --along the cutting plane.
C   --
C   --Parameters:
C   --   CUTPLA - IN - the 3 points defining the cutting plane;
C   --      CUTPLA(1,i) is the x coordinate for point i, etc
C   --   CUTMAT - OUT - the rotation matrix
C   --   * - return statement if the 3 points do not define a plane

      REAL CUTPLA(3,3)
      REAL CUTMAT(3,3)

      X1 = CUTPLA(1,1)
      Y1 = CUTPLA(2,1)
      Z1 = CUTPLA(3,1)
      X2 = CUTPLA(1,2)
      Y2 = CUTPLA(2,2)
      Z2 = CUTPLA(3,2)
      X3 = CUTPLA(1,3)
      Y3 = CUTPLA(2,3)
      Z3 = CUTPLA(3,3)

      X12 = SQRT ((X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)
      X13 = SQRT ((X3-X1)**2 + (Y3-Y1)**2 + (Z3-Z1)**2)
      X23 = SQRT ((X2-X3)**2 + (Y2-Y3)**2 + (Z2-Z3)**2)
      IF (X12 .EQ. 0.0) GOTO 100
      X15 = (X12**2 + X13**2 - X23**2) / (2.0*X12)
      XLAM = (X2-X1) / X12
      YLAM = (Y2-Y1) / X12
      ZLAM = (Z2-Z1) / X12
      X5 = X1 + XLAM*X15
      Y5 = Y1 + YLAM*X15
      Z5 = Z1 + ZLAM*X15
      X35 = SQRT ((X3-X5)**2 + (Y3-Y5)**2 + (Z3-Z5)**2)
      IF (X35 .EQ. 0.0) GOTO 100
      XPSI = (X3-X5) / X35
      YPSI = (Y3-Y5) / X35
      ZPSI = (Z3-Z5) / X35
      XNU = YLAM*ZPSI - ZLAM*YPSI
      YNU = ZLAM*XPSI - XLAM*ZPSI
      ZNU = XLAM*YPSI - YLAM*XPSI

      CUTMAT(1,1) = XLAM
      CUTMAT(2,1) = YLAM
      CUTMAT(3,1) = ZLAM
      CUTMAT(1,2) = XPSI
      CUTMAT(2,2) = YPSI
      CUTMAT(3,2) = ZPSI
      CUTMAT(1,3) = XNU
      CUTMAT(2,3) = YNU
      CUTMAT(3,3) = ZNU

      RETURN

  100 CONTINUE
      CALL PRTERR ('CMDERR', 'Points do not define a plane')
      RETURN 1
      END
