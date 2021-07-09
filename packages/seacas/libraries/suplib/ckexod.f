C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CKEXOD (EXODUS, *)
C=======================================================================

C   --*** CKEXOD *** (ETCLIB) Check for EXODUS format
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --CKEXOD prints an error message if the database is not in the EXODUS
C   --database format.
C   --
C   --Parameters:
C   --   EXODUS - IN - true iff EXODUS format
C   --   * - return statement if error

      LOGICAL EXODUS

      IF (.NOT. EXODUS) THEN
         CALL PRTERR ('CMDERR',
     &      'Command not allowed on GENESIS database')
         RETURN 1
      END IF

      RETURN
      END
