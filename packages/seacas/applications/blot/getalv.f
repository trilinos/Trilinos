C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETALV (A, NALVAR, ALIVAL, ISTEP, LENE, ISEVOK,
     *  ALIVE, VAR)
C=======================================================================

C   --*** GETALV *** (MESH) Read birth/death variable
C   --   Written by Amy Gilkey - revised 10/28/87
C   --
C   --GETALV reads the values for the requested birth/death variable and
C   --returns the element state.
C   --
C   --The element is alive iff the input variable is 0.0.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NALVAR - IN - the variable sequence number
C   --   ALIVAL - IN - the value to indicate element is fully alive
C   --   ISTEP - IN - the time step number
C   --   LENE - IN - the cumulative element counts by element block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   ALIVE - OUT - true iff the element i is alive
C   --   VAR - SCRATCH - the birth/death variable array; may be ALIVE
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      INTEGER LENE(0:*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      LOGICAL ALIVE(NUMEL)
      REAL VAR(NUMEL)

      CHARACTER CDUM

      CALL GETVAR (A, NALVAR, -1, ISTEP, NUMEL, VAR)

      CALL DBVTYP_BL (NALVAR, CDUM, IDALV)
      DO 120 IELB = 1, NELBLK
         IF (ISEVOK(IELB,IDALV)) THEN
            DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
               ALIVE(IEL) = (VAR(IEL) .EQ. ALIVAL)
  100       CONTINUE
         ELSE
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               ALIVE(IEL) = .TRUE.
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
