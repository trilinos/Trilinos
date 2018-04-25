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

C $Log: msscal.f,v $
C Revision 1.5  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.4  2007/11/14 20:14:53  gdsjaar
C Added optional 'alive value' to the death on variable command.  The
C default value is 0.0, but you can now specify a different value to
C indicate aliveness (for example, the presto DEATH_DUMMY_VAR treats 1.0
C as the alive value).
C
C Example: DEATH ON DEATH_DUMMY_VAR 1
C
C Revision 1.3  1997/10/23 12:59:33  gdsjaar
C Fixed initialization order problem that was affecting zooms.
C
C Revision 1.2  1997/09/02 14:52:46  caforsy
C Changed name from HEX_SHELL to HEXSHELL.  Fixed hexshell bugs
C
C Revision 1.1  1994/04/07 20:05:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:04  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE MSSCAL (DOSCAL, NNPSEL, NPSEL, XN, YN, ZN)
C=======================================================================

C   --*** MSSCAL *** (MESH) Set the mesh scale
C   --   Written by Amy Gilkey - revised 08/31/87
C   --
C   --MSSCAL sets the mesh window from the selected nodes, if requested,
C   --and calls SETUP to calculate the mesh scaling variables.
C   --
C   --Parameters:
C   --   DOSCAL - IN - true iff the mesh window limits are to be calculated
C   --   NNPSEL - IN - the number of nodes in the mesh window
C   --   NPSEL - IN - the nodes in the mesh window
C   --   XN, YN, ZN - IN - the nodal coordinates (ZN for 3D only)
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/
C   --   Uses XISSYM, YISSYM, LFTSYM, BOTSYM of /VIEWS/
C   --   Sets and uses ZMMESH of /MSHLIM/
C   --   Uses SQMESH, RDMESH of /MSHLIM/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'mshopt.blk'
      COMMON /VIEWS/  MULTIM,
     &   XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM
      LOGICAL MULTIM, XISSYM, YISSYM, LFTSYM, BOTSYM
      COMMON /MSHLIM/ UNMESH(KFAR), ALMESH(KFAR),
     &   ZMMESH(KTOP), RDMESH(KTOP), TICMSH, SQMESH
      LOGICAL SQMESH
      COMMON /MSHLIC/ MSCTYP
      CHARACTER*8 MSCTYP

      LOGICAL DOSCAL, INHBF
      SAVE INHBF
      INTEGER NPSEL(*)
      REAL XN(*), YN(*), ZN(*)
      DATA INHBF /.TRUE./

C   --Calculate the axis limits

      INHBF = (INHBF .AND.
     *  RDMESH(KLFT) .EQ. 0.0 .and. RDMESH(KRGT) .EQ. 0.0 .AND.
     *  RDMESH(KBOT) .EQ. 0.0 .and. RDMESH(KTOP) .EQ. 0.0)
      IF (DOSCAL .OR. INHBF) THEN
         INHBF = .FALSE.

C      --Calculate the axis limits

         CALL MINMXS (NNPSEL, NPSEL, XN, RDMESH(KLFT), RDMESH(KRGT))
         CALL MINMXS (NNPSEL, NPSEL, YN, RDMESH(KBOT), RDMESH(KTOP))

C      --Expand limits a little

         CALL EXPLIM (2, RDMESH, RDMESH)
         CALL ADJLIM (MSHDEF,
     &      XISSYM, YISSYM, LFTSYM, BOTSYM, XAXSYM, YAXSYM,
     &      SQMESH, RDMESH, ZMMESH)
      END IF

C   --Set up graphics

      CALL SETUP (MSHDEF, ZMMESH)

      RETURN
      END
