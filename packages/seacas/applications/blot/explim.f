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

C $Log: explim.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:00:32  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:54  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE EXPLIM (NDIM, RDMESH, EXMESH)
C=======================================================================

C   --*** EXPLIM *** (BLOT) Expand 2D or 3D mesh limits by 5%
C   --   Written by Amy Gilkey - revised 06/30/86
C   --
C   --EXPLIM expands the mesh limits by 5% of the limits of the maximum
C   --dimension (2.5% on each side).
C   --
C   --Parameters:
C   --   NDIM - IN - the number of dimensions to be expanded
C   --   RDMESH - IN - the mesh limits
C   --      (left, right, bottom, top, near, far)
C   --   EXMESH - OUT - the expanded mesh limits (may be RDMESH)
C   --      (left, right, bottom, top, near, far)

      PARAMETER (PCT2 = 0.025)
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      REAL RDMESH(2*NDIM), EXMESH(2*NDIM)

      DIF = RDMESH(2) - RDMESH(1)
      DO 100 I = 2+2, 2*NDIM, 2
         DIF = MAX (DIF, (RDMESH(I) - RDMESH(I-1)))
  100 CONTINUE
      IF (DIF .EQ. 0.0) THEN
         DO 110 I = 2, 2*NDIM, 2
            DIF = MAX (DIF, ABS (RDMESH(I)))
  110    CONTINUE
         IF (DIF .EQ. 0.0) DIF = 1.0
      END IF
      DIF = PCT2 * DIF

      DO 120 I = 2, 2*NDIM, 2
         EXMESH(I-1) = RDMESH(I-1) - DIF
         EXMESH(I) = RDMESH(I) + DIF
  120 CONTINUE

      RETURN
      END
