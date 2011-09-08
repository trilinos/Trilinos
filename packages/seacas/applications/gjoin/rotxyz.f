C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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

C $Id: rotxyz.f,v 1.1 1999/01/18 19:21:26 gdsjaar Exp $
C $Log: rotxyz.f,v $
C Revision 1.1  1999/01/18 19:21:26  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:28  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1  1992/11/11 22:00:13  gdsjaar
C Added revolve and revcen transformation commands.
C
c Revision 1.1.1.1  1990/11/12  16:10:50  gdsjaar
c GREPOS: Genesis Repositioning
c
c Revision 1.1  90/11/12  16:10:49  gdsjaar
c Initial revision
c 
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
