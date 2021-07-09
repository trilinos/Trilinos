C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSDONE (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS)
C=======================================================================

C   --*** MSDONE *** (MESH) Clean up after mesh plot set
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --MSDONE cleans up the memory after a mesh plot set.  Many variables
C   --may be needed to do "picks" on the mesh plot.
C   --
C   --Parameters:
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

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS

C   --Release the memory reserved for plot set

      IF (IS3DIM .AND. DOIXF) THEN
         CALL MDDEL ('IXFAC')
      END IF
      IF (DODEAD) THEN
         CALL MDDEL ('IDN2B')
      END IF
      IF (DONPS) THEN
         CALL MDDEL ('SCRNPS')
      END IF
      IF (DOESS) THEN
         CALL MDDEL ('SCRESS')
      END IF

      RETURN
      END
