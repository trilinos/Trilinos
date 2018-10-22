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

C $Log: roteye.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:10:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:38  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ROTEYE (EYE, ROTCEN, ROTMAT, *)
C=======================================================================

C   --*** ROTEYE *** (MESH) Get rotation matrix from eye position
C   --   Written by Amy Gilkey - revised 10/21/86
C   --
C   --ROTEYE derives the rotation matrix from the given eye position.
C   --Rotation matrix gives an X rotation followed by a Y rotation.
C   --
C   --Parameters:
C   --   EYE - IN - the eye position
C   --   ROTCEN - IN - the center of rotation
C   --   ROTMAT - OUT - the rotation matrix
C   --   * - return statement if error in eye position

      REAL EYE(3)
      REAL ROTCEN(3)
      REAL ROTMAT(3,3)

      X = EYE(1) - ROTCEN(1)
      Y = EYE(2) - ROTCEN(2)
      Z = EYE(3) - ROTCEN(3)
      VMAG1 = SQRT (Y*Y + Z*Z)
      VMAG2 = SQRT (X*X + Y*Y + Z*Z)
      IF (VMAG1 .EQ. 0.0) GOTO 100
      COS1 = Z / VMAG1
      SIN1 = Y / VMAG1
      COS2 = VMAG1 / VMAG2
      SIN2 = -X / VMAG2

      ROTMAT(1,1) = COS2
      ROTMAT(2,1) = SIN1 * SIN2
      ROTMAT(3,1) = COS1 * SIN2
      ROTMAT(1,2) = 0.0
      ROTMAT(2,2) = COS1
      ROTMAT(3,2) = - SIN1
      ROTMAT(1,3) = - SIN2
      ROTMAT(2,3) = SIN1 * COS2
      ROTMAT(3,3) = COS1 * COS2

      RETURN

  100 CONTINUE
      CALL PRTERR ('CMDERR', 'Eye cannot be exactly on center')
      RETURN 1
      END
