C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSSTEP (A, ISTEP, NEWSET, NEWFAC, ANYDEF, ANYUND,
     &  DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &  LENF, NLNKF, LINKF, LENL, KLNSET,
     &  KHIDEN, KHIDEF, KIXFAC,
     &  KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &  KXF, KYF, KZF, KDXF, KDYF, KDZF,
     &  NEWELB, IELBST, KN2ELB, KDN2B, IF2EL, KNPSUR)
C=======================================================================

C   --*** MSSTEP *** (MESH) Compute coordinates, etc for time step
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --MSSTEP computes the values needed for the mesh plot of a particular
C   --time step.  These values include the nodal deformed coordinates,
C   --the face coordinates (both deformed and undeformed), and
C   --the line set.  The memory for these values is expanded if needed.
C   --The values which do not change between time steps (e.g., the
C   --undeformed coordinates are not recalculated.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   IF2EL2 - the secondary element number of each face
C   --   IE2ELB - the element block for each element
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the time step number
C   --   NEWSET - IN - true iff start of a new plot set
C   --   NEWFAC - IN - true iff the faces have changed
C   --   ANYDEF - IN - true iff any deformed mesh is to be plotted
C   --   ANYUND - IN - true iff any undeformed mesh is to be plotted
C   --   DOIXF - IN - true iff the IXFAC array is needed
C   --   DON2B - IN - true iff the IN2ELB array is needed
C   --   DOELED - IN - true iff the deformed element quarilateral centers
C   --      are needed
C   --   DOELEU - IN - true iff the undeformed element quarilateral centers
C   --      are needed
C   --   DODEAD - IN - true iff dead nodes are needed
C   --   DONPS - IN - true iff node set information is needed
C   --   DOESS - IN - true iff side set information is needed
C   --   LENF - IN/OUT - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN/OUT - the connectivity for all faces
C   --   LENL - IN/OUT - the cumulative line counts by element block
C   --   KLNSET - IN/OUT - the dynamic memory index of the sorted line set
C   --   KHIDEN - OUT - the index of HIDENP - sized only
C   --   KHIDEF - OUT - the index of HIDEF - sized only
C   --   KIXFAC - OUT - the index of IXFAC - sized only
C   --   KXN, KYN, KZN - IN/OUT - the indices of XN, YN, ZN -
C   --      the nodal coordinates (ZN for 3D only)
C   --   KHZ, KVT, KPD - IN/OUT - the indices of HZ, VT, PD -
C   --      the undeformed nodal coordinates
C   --      (rotated for 3D, PD for 3D only)
C   --   KDHZ, KDVT, KDPD - IN/OUT - the indices of DHZ, DVT, DPD -
C   --      the deformed nodal coordinates
C   --      (rotated for 3D, DPD for 3D only)
C   --   KXF, KYF, KZF - IN/OUT - the indices of XF, YF, ZF -
C   --      the element face center coordinates
C   --      (rotated for 3D, ZF for 3D only)
C   --   KDXF, KDYF, KDZF - IN/OUT - the indices of DXF, DYF, DZF -
C   --      the deformed face center coordinates
C   --      (rotated for 3D, DZF for 3D only)
C   --   NEWELB - IN/OUT - the new element blocks flag:
C   --      0 = no new element blocks
C   --      1 = new selected element blocks
C   --      2 = new displayed element blocks (implies new selected blocks)
C   --   IELBST - IN/OUT - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   KN2ELB - IN/OUT - the index of IN2ELB -
C   -       the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   KDN2B - IN/OUT - the index of IDN2B -
C   --      the element block for each dead node; dead if >= 0
C   --   IF2EL - IN - the element number of each face
C   --   KNPSUR - IN/OUT - the index of NPSURF -
C   --      the node numbers of the surface nodes or mesh boundary nodes (2D)
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NNPSUR, NUMNPF, LLNSET of /D3NUMS/
C   --   Uses DEFPRO, DEFOK, DFAC of /DEFORM/
C   --   Uses ROTMAT, ROTCEN of /ROTOPT/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /DEFORM/ DEFPRO, DEFOK, DEFFAC, DDFAC, DFAC,
     &  IXDEF, IYDEF, IZDEF
      LOGICAL DEFPRO, DEFOK
      COMMON /ROTOPT/ NEWROT, ROTMAT(3,3), ROTCEN(3), EYE(3)
      LOGICAL NEWROT

      COMMON /SIZES/  NPSIZ
C      --NPSIZ - the size of the nodal coordinate arrays passed to the
C      --   HIDDEN routine for partial line

      DIMENSION A(*)
      LOGICAL NEWSET, NEWFAC
      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS
      INTEGER LENF(0:NELBLK+4)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER LENL(-2:NELBLK)
      INTEGER NEWELB
      INTEGER IELBST(NELBLK)
      INTEGER IF2EL(*)

      LOGICAL NEWNPD, NEWNPU, NEWELD, NEWELU
      SAVE NEWNPU, NEWELU

      DATA NEWNPU, NEWELU / .TRUE., .TRUE. /

      IF (NEWFAC) THEN
C      --Force the re-calculation of connected element count, if needed
        NEWELB = 1

C      --Force the re-calculation of undeformed nodal and face coordinates
        NEWNPU = .TRUE.
        NEWELU = .TRUE.

C      --Set size of deformed coordinate arrays for hidden line algorithm
        NPSIZ = NUMNPF
        IF (IS3DIM) NPSIZ = (13 * NUMNPF) / 10
      END IF

C   --Force the re-calculation of deformed nodal and face coordinates
      NEWNPD = DEFOK .AND. DEFPRO
      NEWELD = DEFOK .AND. DEFPRO

C   --Force the re-calculation of undeformed nodal and face coordinates
      IF (NEWROT .OR. NEWSET) THEN
        NEWNPU = .TRUE.
        NEWELU = .TRUE.
      END IF

C   --Release memory that will vary in size or not be needed

      IF (NEWSET .OR. NEWFAC) THEN
        IF (NEWELB .GE. 1) THEN
          CALL MDLONG ('IN2ELB', KN2ELB, 0)
        END IF
        IF (DODEAD) CALL MDLONG ('IDN2B', KDN2B, 0)
      END IF
      IF (NNPSUR .LT. 0) THEN
        CALL MDLONG ('NPSURF', KNPSUR, 0)
      END IF
      IF (IS3DIM .AND. NEWNPU) THEN
        CALL MDLONG ('HZ', KHZ, 0)
        CALL MDLONG ('VT', KVT, 0)
        IF (IS3DIM) CALL MDLONG ('PD', KPD, 0)
      END IF
      IF (NEWNPD) THEN
        CALL MDLONG ('DHZ', KDHZ, 0)
        CALL MDLONG ('DVT', KDVT, 0)
        IF (IS3DIM) CALL MDLONG ('DPD', KDPD, 0)
      END IF
      IF (NEWELU) THEN
        CALL MDLONG ('XF', KXF, 0)
        CALL MDLONG ('YF', KYF, 0)
        IF (IS3DIM) CALL MDLONG ('ZF', KZF, 0)
      END IF
      IF (NEWELD) THEN
        CALL MDLONG ('DXF', KDXF, 0)
        CALL MDLONG ('DYF', KDYF, 0)
        IF (IS3DIM) CALL MDLONG ('DZF', KDZF, 0)
      END IF
      IF (IS3DIM .AND. NEWFAC) THEN
        CALL MDLONG ('HIDENP', KHIDEN, 0)
        CALL MDLONG ('HIDEF', KHIDEF, 0)
        IF (DOIXF) CALL MDLONG ('IXFAC', KIXFAC, 0)
      END IF

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Compute new line set if new faces

      IF (NEWFAC) THEN
        CALL MSLINS (A, LENF, NLNKF, LINKF, IF2EL, LENL, KLNSET)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 100
      END IF

C   --Compute node to element block and dead node to element block,
C   --if needed

      IF ((NEWELB .GE. 1) .AND. (DON2B .OR. DONPS .OR. DOESS)) THEN
        IF (NEWSET .OR. NEWFAC) THEN
          CALL MDLONG ('IN2ELB', KN2ELB, NUMNPF)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100
        END IF

        CALL MAKN2B (LENF, NLNKF, LINKF, IELBST, A(KN2ELB))
        NEWELB = 0
      END IF

      IF (DODEAD) THEN
        IF (NEWSET .OR. NEWFAC) THEN
          CALL MDLONG ('IDN2B', KDN2B, NUMNPF)
        END IF
        IF (IS3DIM) THEN
          CALL MDFIND ('IF2EL2', KIF2E2, IDUM)
        ELSE
          KIF2E2 = 1
          KE2ELB = 1
        END IF
        CALL MDFIND ('IE2ELB', KE2ELB, IDUM)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 100

        CALL MAKD2B (LENF, NLNKF, LINKF, IELBST,
     &    IF2EL, A(KIF2E2), A(KE2ELB), A(KDN2B))

      END IF

C   --Compute NPSURF nodes for 3D, if changed

      IF (NNPSUR .LT. 0) THEN
        CALL MDLONG ('NPSURF', KNPSUR, NUMNPF)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 100

        CALL MAKSUR (LENF, NLNKF, LINKF,
     &    DODEAD, A(KDN2B), A(KNPSUR))

        CALL MDLONG ('NPSURF', KNPSUR, NNPSUR)
      END IF

C   --Put rotated undeformed coordinates in undeformed arrays, if needed

      IF (NEWNPU .AND. ANYUND) THEN
        IF (IS3DIM) THEN
          CALL MDLONG ('HZ', KHZ, NPSIZ)
          CALL MDLONG ('VT', KVT, NPSIZ)
          IF (IS3DIM) CALL MDLONG ('PD', KPD, NPSIZ)
        ELSE
          KHZ = KXN
          KVT = KYN
          KPD = KZN
        END IF

        IF (IS3DIM) THEN
          CALL BL_ROTATE (NNPSUR, A(KNPSUR), ROTMAT, ROTCEN,
     &      A(KXN), A(KYN), A(KZN), A(KHZ), A(KVT), A(KPD))

C         --Force the re-calculation of undeformed face coordinates
          NEWELU = .TRUE.
        END IF

        NEWNPU = .FALSE.
        NEWROT = .FALSE.
      END IF

C   --Compute deformed nodal coordinates, if needed

      IF (NEWNPD .AND. ANYDEF) THEN
        CALL MDLONG ('DHZ', KDHZ, NPSIZ)
        CALL MDLONG ('DVT', KDVT, NPSIZ)
        IF (IS3DIM) CALL MDLONG ('DPD', KDPD, NPSIZ)

        CALL DEFXYZ (A, ISTEP, DFAC, IS3DIM, A(KNPSUR),
     &    A(KXN), A(KYN), A(KZN), A(KDHZ), A(KDVT), A(KDPD))

        IF (IS3DIM) THEN
          CALL BL_ROTATE (NNPSUR, A(KNPSUR), ROTMAT, ROTCEN,
     &      A(KDHZ), A(KDVT), A(KDPD), A(KDHZ), A(KDVT), A(KDPD))
        END IF

C      --Force the re-calculation of deformed face coordinates
        NEWELD = .TRUE.

        NEWNPU = .FALSE.
      END IF

C   --Get face coordinates, if changed

      IF (NEWELU .AND. DOELEU) THEN
        IF (NEWELU) THEN
          CALL MDLONG ('XF', KXF, LENF(NELBLK))
          CALL MDLONG ('YF', KYF, LENF(NELBLK))
          IF (IS3DIM) CALL MDLONG ('ZF', KZF, LENF(NELBLK))
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100
        END IF

        CALL ELECOR (NDIM, NELBLK, LENF, NLNKF, LINKF,
     &    A(KHZ), A(KVT), A(KPD), A(KXF), A(KYF), A(KZF))
        NEWELU = .FALSE.
      END IF

      IF (NEWELD .AND. DOELED) THEN
        IF (NEWELD) THEN
          CALL MDLONG ('DXF', KDXF, LENF(NELBLK))
          CALL MDLONG ('DYF', KDYF, LENF(NELBLK))
          IF (IS3DIM) CALL MDLONG ('DZF', KDZF, LENF(NELBLK))
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100
        END IF

        CALL ELECOR (NDIM, NELBLK, LENF, NLNKF, LINKF,
     &    A(KDHZ), A(KDVT), A(KDPD), A(KDXF), A(KDYF), A(KDZF))
        NEWELD = .FALSE.
      END IF

      IF (IS3DIM .AND. NEWFAC) THEN
        CALL MDLONG ('HIDENP', KHIDEN, NUMNPF)
        CALL MDLONG ('HIDEF', KHIDEF, LENF(NELBLK))
        IF (DOIXF) CALL MDLONG ('IXFAC', KIXFAC, LENF(NELBLK))
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 100
      END IF

 100  CONTINUE
      RETURN
      END
