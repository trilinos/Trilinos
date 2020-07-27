C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE QNPICK (TYPE, ISUP, ISDEF,
     &  A, KXN, KYN, KZN, KHIDEN, KNPSUR)
C=======================================================================

C   --*** QNPICK *** (MESH) Query last mesh display for node information
C   --   Written by Amy Gilkey - revised 04/06/88
C   --
C   --QNPICK returns the address of the last displayed mesh characteristics
C   --for cursor input.  This routine assumes that the characteristics
C   --are in dynamic memory.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   XN, YN, ZN - the nodal coordinates (ZN for 3D only)
C   --   DHZ, DVT, DPD - the deformed nodal coordinates (DPD for 3D only)
C   --   HIDENP - true iff node i is hidden (3D only)
C   --   NPSURF - the node numbers of the surface nodes (3D only)
C   --
C   --Parameters:
C   --   TYPE - IN - the type of mesh information needed:
C   --      'QUERY'     - return ISUP and ISDEF only
C   --      'ORIGINAL'  - original mesh coordinates
C   --      'DISPLAYED' - last displayed mesh coordinates
C   --   ISUP - OUT - true iff a mesh is displayed
C   --   ISDEF - OUT - true iff the displayed mesh is deformed
C   --   A - IN - the dynamic memory base array
C   --   KXN, KYN, KZN - OUT - the dynamic memory index of the coordinates
C   --      requested (KZN in 3D only)
C   --   KHIDEN - OUT - the dynamic memory index of the node hidden indicator
C   --      (3D only)
C   --   KNPSUR - OUT - the dynamic memory index of the surface node indices
C   --      (3D only)
C   --
C   --Common Variables:
C   --   Uses /PICK/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      COMMON /PICK/   INITP, PKDEF,
     &  PKMESH(KTOP), PKRMAT(3,3), PKRCEN(3),
     &  DMESH(KTOP), DVRAT, DXMID, DYMID, DXLAST, DYLAST
      LOGICAL INITP, PKDEF

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ISUP, ISDEF
      CHARACTER*(*) TYPE

      ISUP = INITP
      IF (ISUP) THEN
        ISDEF = PKDEF
      ELSE
        ISDEF = .FALSE.
      END IF
      KXN = 1
      KYN = 1
      KZN = 1
      KHIDEN = 1
      KNPSUR = 1

      IF (TYPE .EQ. 'ORIGINAL') THEN
        CALL MDFIND ('XN', KXN, IDUM)
        CALL MDFIND ('YN', KYN, IDUM)
        IF (IS3DIM) CALL MDFIND ('ZN', KZN, IDUM)
      END IF

      IF (ISUP) THEN
        IF (TYPE .EQ. 'DISPLAYED') THEN
          IF (.NOT. PKDEF) THEN
            IF (.NOT. IS3DIM) THEN
              CALL MDFIND ('XN', KXN, IDUM)
              CALL MDFIND ('YN', KYN, IDUM)
              IF (IS3DIM) CALL MDFIND ('ZN', KZN, IDUM)
            ELSE
              CALL MDFIND ('HZ', KXN, IDUM)
              CALL MDFIND ('VT', KYN, IDUM)
              IF (IS3DIM) CALL MDFIND ('PD', KZN, IDUM)
            END IF
          ELSE
            CALL MDFIND ('DHZ', KXN, IDUM)
            CALL MDFIND ('DVT', KYN, IDUM)
            IF (IS3DIM) CALL MDFIND ('DPD', KZN, IDUM)
          END IF
          IF (IS3DIM) CALL MDFIND ('HIDENP', KHIDEN, IDUM)
          IF (IS3DIM) CALL MDFIND ('NPSURF', KNPSUR, IDUM)
        END IF
      END IF

      RETURN
      END
