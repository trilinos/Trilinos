C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C $Log: scal3d.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:10:50  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:57  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCAL3D (MSCTYP, ROTMAT, ROTCEN, ALMESH, D2MESH)
C=======================================================================

C   --*** SCAL3D *** (MESH) Find the 2D limits of a 3D mesh
C   --   Written by Amy Gilkey - revised 06/30/86
C   --
C   --SCAL3D takes the mesh limits and scales them in one of two ways:
C   --   1) it finds the 2D limits enclosing the smallest sphere enclosing
C   --      a cube from the maximum 3D limits
C   --   2) it finds the 2D limits of the rotated mesh
C   --The mesh limits are expanded a little.
C   --
C   --Parameters:
C   --   MSCTYP - IN - mesh scaling flag (as in /MSHLIM/)
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the mesh center for the rotation
C   --   ALMESH - IN - the mesh limits
C   --   D2MESH - OUT - the 2D rotated mesh limits

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      CHARACTER*(*) MSCTYP
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL ALMESH(KFAR), D2MESH(KTOP)

      IF (MSCTYP .EQ. 'ALL') THEN
         RMAX = MAX (ALMESH(KRGT)-ROTCEN(1), ROTCEN(1)-ALMESH(KLFT)
     &      ,        ALMESH(KTOP)-ROTCEN(2), ROTCEN(2)-ALMESH(KBOT)
     &      ,        ALMESH(KFAR)-ROTCEN(3), ROTCEN(3)-ALMESH(KNEA))
         RAD = SQRT(3.0) * RMAX
         D2MESH(KLFT) = - RAD
         D2MESH(KRGT) = + RAD
         D2MESH(KBOT) = - RAD
         D2MESH(KTOP) = + RAD

      ELSE
         CALL BL_ROTATE (1, 1, ROTMAT, ROTCEN,
     &      ALMESH(KLFT), ALMESH(KBOT), ALMESH(KNEA), X1, Y1, RDUM)
         CALL BL_ROTATE (1, 1, ROTMAT, ROTCEN,
     &      ALMESH(KRGT), ALMESH(KTOP), ALMESH(KFAR), X2, Y2, RDUM)
         D2MESH(KLFT) = MIN (X1, X2)
         D2MESH(KRGT) = MAX (X1, X2)
         D2MESH(KBOT) = MIN (Y1, Y2)
         D2MESH(KTOP) = MAX (Y1, Y2)
      END IF

      CALL EXPLIM (2, D2MESH, D2MESH)

      RETURN
      END
