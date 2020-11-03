C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INPICK (TYPE)
C=======================================================================

C   --*** INPICK *** (MESH) Initialize cursor input for last display
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --INPICK saves the displayed mesh characteristics for cursor input
C   --in common /PICK/.  This routine should be called after the last
C   --mesh display before any cursor input.
C   --
C   --Parameters:
C   --   TYPE - IN - the type of mesh last displayed (determines which
C   --      coordinates arrays to use):
C   --      'DEFORM'   - deformed mesh (or any 3D mesh)
C   --      'UNDEFORM' - undeformed mesh
C   --      'NONE'     - no mesh displayed
C   --
C   --Common Variables:
C   --   Sets /PICK/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses ZMMESH of /MSHLIM/
C   --   Uses ROTMAT, ROTCEN of /ROTOPT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      COMMON /PICK/   INITP, PKDEF,
     &   PKMESH(KTOP), PKRMAT(3,3), PKRCEN(3),
     &   DMESH(KTOP), DVRAT, DXMID, DYMID, DXLAST, DYLAST
      LOGICAL INITP, PKDEF

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /MSHLIM/ UNMESH(KFAR), ALMESH(KFAR),
     &   ZMMESH(KTOP), RDMESH(KTOP), TICMSH, SQMESH
      LOGICAL SQMESH
      COMMON /ROTOPT/ NEWROT, ROTMAT(3,3), ROTCEN(3), EYE(3)
      LOGICAL NEWROT

      CHARACTER*(*) TYPE

      INITP = .FALSE.

      IF ((TYPE .EQ. 'DEFORM') .OR. (TYPE .EQ. 'UNDEFORM')) THEN

C      --Set the pick initialization flag
         INITP = .TRUE.
         PKDEF = (TYPE .EQ. 'DEFORM')

C      --Save the limits of the displayed mesh
         CALL CPYREA (KTOP, ZMMESH, PKMESH)

C      --Save the rotation matrix and the rotation center of the displayed mesh
         IF (IS3DIM) THEN
            CALL CPYREA (3*3, ROTMAT, PKRMAT)
            CALL CPYREA (3, ROTCEN, PKRCEN)
         END IF

C      --Calculate the device coordinates of the limits of the displayed mesh
         CALL MP2PT (1, PKMESH(KLFT), PKMESH(KBOT),
     &      DMESH(KLFT), DMESH(KBOT), M)
         CALL MP2PT (1, PKMESH(KRGT), PKMESH(KTOP),
     &      DMESH(KRGT), DMESH(KTOP), M)

C      --Calculate the device to mesh conversion factor
         DVRAT = (PKMESH(KRGT) - PKMESH(KLFT))
     &      / (DMESH(KRGT) - DMESH(KLFT))

C      --Find the device coordinates of the center of the displayed mesh,
C      --for use as the default starting point
         PKXMID = 0.5 * (PKMESH(KRGT) + PKMESH(KLFT))
         PKYMID = 0.5 * (PKMESH(KTOP) + PKMESH(KBOT))
         CALL MP2PT (1, PKXMID, PKYMID, DXMID, DYMID, M)

         DXLAST = DXMID
         DYLAST = DYMID

      ELSE

C      --Reset the pick initialization flag
         INITP = .FALSE.
      END IF

      RETURN
      END
