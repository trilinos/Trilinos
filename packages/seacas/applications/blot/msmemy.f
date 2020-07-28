C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSMEMY (A, ANYDEF, ANYUND,
     &  DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &  NSURFS, KN2ELB, KDN2B, KIF2EL,
     &  KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &  KXF, KYF, KZF, KDXF, KDYF, KDZF, KHIDEN, KHIDEF, KIXFAC,
     &  KSNPS, KSESS)
C=======================================================================

C   --*** MSMEMY *** (MESH) Reserve memory for mesh plot
C   --   Written by Amy Gilkey - revised 07/27/88
C   --
C   --MSMEMY reserves or finds the memory needed to plot a mesh.
C   --
C   --This routine reserves and finds the following dynamic memory arrays:
C   --   IN2ELB - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   IDN2B - the element block for each dead node; dead if >= 0
C   --   HZ, VT, PD - the undeformed nodal coordinates
C   --      (rotated for 3D, PD for 3D only)
C   --   DHZ, DVT, DPD - the deformed nodal coordinates
C   --      (rotated for 3D, DPD for 3D only)
C   --   XF, YF, ZF - the undeformed face center coordinates
C   --      (rotated for 3D, ZF for 3D only)
C   --   DXF, DYF, DZF - the deformed face center coordinates
C   --      (rotated for 3D, DZF for 3D only)
C   --   HIDENP(i) - true iff node i is hidden (3D only)
C   --   HIDEF(i) - true iff face i is hidden (3D only)
C   --   IXFAC - the indices of the ordered faces (3D only)
C   --   IX2NP - the node number for each mesh index (MASTER process only)
C   --   SCRNPS - scratch array for node sets
C   --   SCRESS - scratch array for side sets
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   IF2EL - the element number of each face
C   --   XN, YN, ZN - IN - the nodal coordinates (ZN for 3D only)
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
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
C   --   NSURFS - IN - the number of surface faces
C   --   KN2ELB - OUT - the index of IN2ELB
C   --   KDN2B - OUT - the index of IDN2B
C   --   KIF2EL - OUT - the index of IF2EL
C   --   KXN, KYN, KZN - OUT - the indices of XN, YN, ZN
C   --   KHZ, KVT, KPD - OUT - the indices of HZ, VT, PD
C   --   KDHZ, KDVT, KDPD - OUT - the indices of DHZ, DVT, DPD
C   --   KXF, KYF, KZF - OUT - the indices of XF, YF, ZF
C   --   KDXF, KDYF, KDZF - OUT - the indices of DXF, DYF, DZF
C   --   KHIDEN - OUT - the index of HIDENP
C   --   KHIDEF - OUT - the index of HIDEF
C   --   KIXFAC - OUT - the index of IXFAC
C   --   KSNPS - OUT - the index of SCRNPS
C   --   KSESS - OUT - the index of SCRESS
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses DEFPRO, DEFOK of /DEFORM/
C   --   Sets NPSIZ of /SIZES/

      include 'dbnumgq.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'sizes.blk'

      DIMENSION A(*)
      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
        FIRST = .FALSE.

C      --Reserve memory (most will be expanded when needed)

        CALL MDRSRV ('IN2ELB', KN2ELB, 0)
        CALL MDFIND ('IF2EL', KIF2EL, IDUM)

        CALL MDFIND ('XN', KXN, IDUM)
        CALL MDFIND ('YN', KYN, IDUM)
        IF (IS3DIM) THEN
          CALL MDFIND ('ZN', KZN, IDUM)
        ELSE
          KZN = 1
        END IF

C      --Set size of deformed coordinate arrays for hidden line algorithm
        NPSIZ = NUMNPF
        IF (IS3DIM) NPSIZ = (13 * NUMNPF) / 10

        IF (IS3DIM) THEN
C         --QNPICK uses MDFIND to find HZ, VT, PD
          CALL MDRSRV ('HZ', KHZ, 0)
          CALL MDRSRV ('VT', KVT, 0)
          CALL MDRSRV ('PD', KPD, 0)
        ELSE
          KHZ = KXN
          KVT = KYN
          KPD = KZN
        END IF

        IF (DEFOK) THEN
C         --QNPICK uses MDFIND to find DHZ, DVT, DPD
          CALL MDRSRV ('DHZ', KDHZ, 0)
          CALL MDRSRV ('DVT', KDVT, 0)
          IF (IS3DIM) THEN
            CALL MDRSRV ('DPD', KDPD, 0)
          ELSE
            KDPD = 1
          END IF
        ELSE
          KDHZ = 1
          KDVT = 1
          KDPD = 1
        END IF

        CALL MDRSRV ('XF', KXF, 0)
        CALL MDRSRV ('YF', KYF, 0)
        IF (IS3DIM) THEN
          CALL MDRSRV ('ZF', KZF, 0)
        ELSE
          KZF = 1
        END IF

        IF (DEFOK) THEN
          CALL MDRSRV ('DXF', KDXF, 0)
          CALL MDRSRV ('DYF', KDYF, 0)
          IF (IS3DIM) THEN
            CALL MDRSRV ('DZF', KDZF, 0)
          ELSE
            KDZF = 1
          END IF
        END IF

C      --QNPICK uses MDFIND to find HIDENP
        IF (IS3DIM) THEN
          CALL MDRSRV ('HIDENP', KHIDEN, NUMNPF)
        ELSE
          KHIDEN = 1
        END IF
        IF (IS3DIM) THEN
          CALL MDRSRV ('HIDEF', KHIDEF, NSURFS)
        ELSE
          KHIDEF = 1
        END IF

C      --Reserve the mesh index to node number index
        CALL MDRSRV ('IX2NP', KIX2NP, 0)
      END IF

C   --Find memory

      CALL MDFIND ('IN2ELB', KN2ELB, IDUM)
      CALL MDFIND ('IF2EL', KIF2EL, IDUM)

      IF (DODEAD) THEN
        CALL MDRSRV ('IDN2B', KDN2B, 0)
      ELSE
        KDN2B = 1
      END IF

      CALL MDFIND ('XN', KXN, IDUM)
      CALL MDFIND ('YN', KYN, IDUM)
      IF (IS3DIM) THEN
        CALL MDFIND ('ZN', KZN, IDUM)
      ELSE
        KZN = 1
      END IF

      IF (IS3DIM) THEN
        CALL MDFIND ('HZ', KHZ, IDUM)
        CALL MDFIND ('VT', KVT, IDUM)
        CALL MDFIND ('PD', KPD, IDUM)
      ELSE
        KHZ = KXN
        KVT = KYN
        KPD = KZN
      END IF

      IF (DEFOK .AND. DEFPRO) THEN
        CALL MDFIND ('DHZ', KDHZ, IDUM)
        CALL MDFIND ('DVT', KDVT, IDUM)
        IF (IS3DIM) THEN
          CALL MDFIND ('DPD', KDPD, IDUM)
        ELSE
          KDPD = 1
        END IF
      ELSE
        KDHZ = 1
        KDVT = 1
        KDPD = 1
      END IF

      CALL MDFIND ('XF', KXF, IDUM)
      CALL MDFIND ('YF', KYF, IDUM)
      IF (IS3DIM) THEN
        CALL MDFIND ('ZF', KZF, IDUM)
      ELSE
        KZF = 1
      END IF

      IF (DEFOK .AND. DEFPRO) THEN
        CALL MDFIND ('DXF', KDXF, IDUM)
        CALL MDFIND ('DYF', KDYF, IDUM)
        IF (IS3DIM) THEN
          CALL MDFIND ('DZF', KDZF, IDUM)
        ELSE
          KDZF = 1
        END IF
      ELSE
        KDXF = 1
        KDYF = 1
        KDZF = 1
      END IF

      IF (IS3DIM) THEN
        CALL MDFIND ('HIDENP', KHIDEN, IDUM)
      ELSE
        KHIDEN = 1
      END IF
      KIXFAC = 1
      IF (IS3DIM) THEN
        CALL MDFIND ('HIDEF', KHIDEF, IDUM)
        IF (DOIXF) CALL MDRSRV ('IXFAC', KIXFAC, NSURFS)
      ELSE
        KHIDEF = 1
      END IF

      IF (DONPS) THEN
        CALL MDRSRV ('SCRNPS', KSNPS, NUMNPS)
      ELSE
        KSNPS = 1
      END IF
      IF (DOESS) THEN
        CALL MDRSRV ('SCRESS', KSESS, NUMESS)
      ELSE
        KSESS = 1
      END IF

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

 100  CONTINUE
      RETURN
      END
