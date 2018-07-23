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

C $Log: scalax.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:10:54  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:59  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCALAX
C=======================================================================

C   --*** SCALAX *** (MESH) Set axis scale
C   --   Written by Amy Gilkey - revised 01/29/88
C   --
C   --SCALAX sets the new zoom window limits based on the scale type.
C   --
C   --Parameters:
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DFAC of /DEFORM/
C   --   Uses MSHDEF of /MSHOPT/
C   --   Uses XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM of /VIEWS/
C   --   Uses UNMESH, ALMESH, MSCTYP, SQMESH of /MSHLIM/
C   --   Sets ZMMESH, RDMESH of /MSHLIM/
C   --   Uses ROTMAT, ROTCEN of /ROTOPT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      
      IF ((MSCTYP .NE. 'EACH') .AND. (MSCTYP .NE. 'SELECTED')) THEN

         IF (MSCTYP .EQ. 'ZOOM') THEN
            CONTINUE

         ELSE IF (MSCTYP .EQ. 'MESH') THEN
            IF (DFAC .EQ. 0.0) THEN
               CALL EXPLIM (2, UNMESH, RDMESH)
            ELSE
               CALL EXPLIM (2, ALMESH, RDMESH)
            END IF

         ELSE IF ((MSCTYP .EQ. 'ROTATION')
     &      .OR. (MSCTYP .EQ. 'ALL')) THEN
            IF (DFAC .EQ. 0.0) THEN
               CALL SCAL3D (MSCTYP, ROTMAT, ROTCEN, UNMESH, RDMESH)
            ELSE
               CALL SCAL3D (MSCTYP, ROTMAT, ROTCEN, ALMESH, RDMESH)
            END IF
         END IF

         CALL ADJLIM (MSHDEF,
     &      XISSYM, YISSYM, LFTSYM, BOTSYM, XAXSYM, YAXSYM,
     &      SQMESH, RDMESH, ZMMESH)
      END IF

      RETURN
      END
