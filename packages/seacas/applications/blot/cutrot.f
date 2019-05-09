C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Log: cutrot.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:57:41  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:09  gdsjaar
c Added RCS Id and Log to all files
c
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
