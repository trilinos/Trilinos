C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION CPUIFC (LDUM)
C=======================================================================

C   --*** CPUIFC *** Dummy cancel function
C   --   Written by Amy Gilkey - revised 02/11/88
C   --
C   --CPUIFC returns the cancel flag as false.

      LOGICAL LDUM

      CPUIFC = .FALSE.

      RETURN
      END
