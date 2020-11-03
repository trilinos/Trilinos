C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPTRND (A, MXSTEP, TYP, NWRDS, NVAR, DATA)
C=======================================================================

C   --*** SPTRND *** (SPLOT) Transfer variables to random file
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --SPTRND transfers needed variables from the sequential database to the
C   --random file (using GETVAR).
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   MXSTEP - IN - the maximum time step number needed
C   --   TYP - IN - the variable type ("N"odal or "E"lement)
C   --   NWRDS - IN - the number of words in a data record
C   --   NVAR - IN - the number of records to read
C   --   DATA - SCRATCH - size = NWRDS
C   --
C   --Common Variables:
C   --   Uses NSPVAR, ISVID of /SPVARS/

      include 'params.blk'
      include 'spvars.blk'

      DIMENSION A(*)
      CHARACTER TYP
      REAL DATA(NWRDS)

      LOGICAL NEED

      CALL DBVIX_BL (TYP, 1, ISID)
      CALL DBVIX_BL (TYP, NVAR, IEID)

      DO 100 ID = ISID, IEID

C      --Determine if variable is needed

         NEED = (LOCINT (ID, NSPVAR, ISVID) .GT. 0)

C      --If so, transfer variable to random file

         IF (NEED) THEN
            CALL GETVAR (A, ID, -1, -MXSTEP, NWRDS, DATA)
         END IF
  100 CONTINUE

      RETURN
      END
