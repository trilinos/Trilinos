C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EVAROK (NVARS, NVAR, NELBLK, IELBST, ISEVOK, ISVOK)
C=======================================================================

C   --*** EVAROK *** (MESH) Get the multi-variable truth table
C   --   Written by Amy Gilkey - revised 10/29/87
C   --
C   --EVAROK creates the multi-variable truth table.  It uses the element
C   --variable truth table and the selected variables to create a table
C   --for only the selected variables.  If no element variables are given,
C   --the truth table is all true.
C   --
C   --Parameters:
C   --   NVARS - IN - the number of variables
C   --   NVAR - IN - the variable numbers, if any
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   ISVOK - OUT - the variable truth table; true iff all variables
C   --      of block j exist and are selected

      INTEGER NVAR(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      LOGICAL ISVOK(NELBLK)

      CHARACTER TYP

      DO 100 IELB = 1, NELBLK
         ISVOK(IELB) = (IELBST(IELB) .GT. 0)
  100 CONTINUE

      DO 120 IVAR = 1, NVARS
         CALL DBVTYP_BL (NVAR(IVAR), TYP, ID)
         IF (TYP .EQ. 'E') THEN
            DO 110 IELB = 1, NELBLK
               IF (.NOT. ISEVOK(IELB,ID)) THEN
                  ISVOK(IELB) = .FALSE.
               END IF
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
