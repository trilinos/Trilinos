C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
