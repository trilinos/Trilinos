C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNKFNC (NUMSTO, *)
C=======================================================================

C   --*** LNKFNC *** (ALGEBRA) Assign storage for time functions
C   --   Written by Amy Gilkey - revised 07/22/87
C   --
C   --LNKFNC sets up the storage locations for the time functions that
C   --need storage for results that must be saved over time steps.
C   --
C   --Parameters:
C   --   NUMSTO - IN/OUT - the number of variable storage locations needed
C   --   * - return statement if an error is found; message is printed
C   --
C   --Common Variables:
C   --   Sets ITMENT of /ENT../
C   --   Uses NUMEQN, NUMENT, TYPENT, INXENT of /ENT../
C   --   Uses FNCSTO of /FNCTB./

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_ent.blk'
      include 'ag_fnctbc.blk'

C   --Allocate storage for time functions

      DO 110 NEQN = 1, NUMEQN
         DO 100 NENT = 3, NUMENT(NEQN)
            IF (TYPENT(NENT,NEQN) .EQ. 'F') THEN
               INX = INXENT(NENT,NEQN)
               IF (FNCSTO(INX)) THEN
                  NUMSTO = NUMSTO + 1
                  ITMENT(NENT,NEQN) = NUMSTO
               ELSE
                  ITMENT(NENT,NEQN) = 0
               END IF
            END IF
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
