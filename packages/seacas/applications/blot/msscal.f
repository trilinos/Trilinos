C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
