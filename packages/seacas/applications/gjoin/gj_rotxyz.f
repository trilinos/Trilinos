C Copyright (c) 2008-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C
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
C

C=======================================================================
      SUBROUTINE ROTXYZ (XYZ, ANG, ROTMAT)
C=======================================================================

C   --*** ROTXYZ *** (GEN3D) Multiply rotation matrix by XYZ rotation
C   --   Written by Amy Gilkey - revised 05/23/86
C   --
C   --ROTXYZ multiplies the current rotation matrix by the given X or Y
C   --or Z axis rotation.  Thus the XYZ axes are in viewing space, not
C   --object space.
C   --
C   --Parameters:
C   --   XYZ - IN - the axis of rotation (X,Y,Z)
C   --   ANG - IN - the angle of rotation (in radians)
C   --   ROTMAT - IN/OUT - the rotation matrix

      CHARACTER XYZ
      REAL ROTMAT(3,3)

      REAL BY(3,3), RES(3,3)

      IF (XYZ .EQ. 'X') THEN
         N1 = 2
         N2 = 3
         N3 = 1
      ELSE IF (XYZ .EQ. 'Y') THEN
         N1 = 3
         N2 = 1
         N3 = 2
      ELSE IF (XYZ .EQ. 'Z') THEN
         N1 = 1
         N2 = 2
         N3 = 3
      ELSE
         n1 = 0
         n2 = 0
         n3 = 0
         CALL PRTERR ('ERROR', 'Invalid axis specification in rotxyz')
         return
      END IF

      COSANG = COS (ANG)
      SINANG = SIN (ANG)
      BY(N1,N1) = COSANG
      BY(N2,N1) = -SINANG
      BY(N1,N3) = 0.0
      BY(N1,N2) = SINANG
      BY(N2,N2) = COSANG
      BY(N2,N3) = 0.0
      BY(N3,N1) = 0.0
      BY(N3,N2) = 0.0
      BY(N3,N3) = 1.0

      DO 20 I = 1, 3
         DO 10 J = 1, 3
            RES(I,J) = ROTMAT(I,1)*BY(1,J) + ROTMAT(I,2)*BY(2,J)
     &         + ROTMAT(I,3)*BY(3,J)
   10    CONTINUE
   20 CONTINUE

      CALL CPYREA (3*3, RES, ROTMAT)

      RETURN
      END
