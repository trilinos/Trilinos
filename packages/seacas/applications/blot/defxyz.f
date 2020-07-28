C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DEFXYZ (A, ISTEP, FACTOR, SURONL, NPSURF,
     &   XN, YN, ZN, DXN, DYN, DZN)
C=======================================================================

C   --*** DEFXYZ *** (MESH) Calculate deformed mesh coordinates
C   --   Written by Amy Gilkey - revised 10/12/87
C   --   D. P. Flanagan, 11/17/82
C   --
C   --DEFXYZ calculates the deformed mesh coordinates for a time step.
C   --The displacement variables for the time steps are input and the
C   --deformed mesh coordinates are calculated.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the time step number
C   --   FACTOR - IN - the magnification factor
C   --   SURONL - IN - deform NPSURF nodes only if true (always true for 3D)
C   --   XN, YN, ZN - IN - the nodal coordinates (ZN for 3D only)
C   --   NPSURF - IN - the node numbers of the surface nodes
C   --      or mesh boundary nodes (2D)
C   --   DXN, DYN, DZN - OUT - the deformed nodal coordinates
C   --      (DZN for 3D only)
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Sets IXDEF, IYDEF, IZDEF of /DEFORM/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /DEFORM/ DEFPRO, DEFOK, DEFFAC, DDFAC, DFAC,
     &   IXDEF, IYDEF, IZDEF
      LOGICAL DEFPRO, DEFOK

      DIMENSION A(*)
      LOGICAL SURONL
      INTEGER NPSURF(*)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      REAL DXN(NUMNPF), DYN(NUMNPF), DZN(NUMNPF)

      CALL GTMVAR (A, IXDEF, -999, ISTEP, NUMNPF, DXN)
      CALL GTMVAR (A, IYDEF, -999, ISTEP, NUMNPF, DYN)
      IF (IS3DIM) CALL GTMVAR (A, IZDEF, -999, ISTEP, NUMNPF, DZN)

      IF (.NOT. IS3DIM) THEN
         IF (SURONL) THEN
            DO 100 IX = 1, NNPSUR
               INP = NPSURF(IX)
               DXN(INP) = XN(INP) + FACTOR * DXN(INP)
               DYN(INP) = YN(INP) + FACTOR * DYN(INP)
  100       CONTINUE
         ELSE
            DO 110 INP = 1, NUMNPF
               DXN(INP) = XN(INP) + FACTOR * DXN(INP)
               DYN(INP) = YN(INP) + FACTOR * DYN(INP)
  110       CONTINUE
         END IF

      ELSE
         DO 120 IX = 1, NNPSUR
            INP = NPSURF(IX)
            DXN(INP) = XN(INP) + FACTOR * DXN(INP)
            DYN(INP) = YN(INP) + FACTOR * DYN(INP)
            DZN(INP) = ZN(INP) + FACTOR * DZN(INP)
  120    CONTINUE
      END IF

      RETURN
      END
