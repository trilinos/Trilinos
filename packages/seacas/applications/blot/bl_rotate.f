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

C=======================================================================
      SUBROUTINE BL_ROTATE (NUM, NPROT, ROTMAT, ROTCEN,
     &   XN, YN, ZN, HZ, VT, PD)
C=======================================================================

C   --*** ROTATE *** (MESH) Rotate 3D coordinates
C   --   Written by Amy Gilkey - revised 09/09/87
C   --
C   --ROTATE rotates the 3D coordinates by subtracting the rotation center
C   --and multipling by the rotation matrix.
C   --
C   --Parameters:
C   --   NUM - IN - the number of nodes to rotate
C   --   NPROT - IN - the node numbers of the nodes to rotate
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the center of the rotation
C   --   XN, YN, ZN - IN - the original nodal coordinates
C   --   HZ, VT, PD - OUT - the rotated nodal coordinates

      INTEGER NPROT(NUM)
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL XN(*), YN(*), ZN(*)
      REAL HZ(*), VT(*), PD(*)

      DO 100 IX = 1, NUM
         INP = NPROT(IX)
         X = XN(INP) - ROTCEN(1)
         Y = YN(INP) - ROTCEN(2)
         Z = ZN(INP) - ROTCEN(3)
         HZ(INP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
         VT(INP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
         PD(INP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
  100 CONTINUE

      RETURN
      END
