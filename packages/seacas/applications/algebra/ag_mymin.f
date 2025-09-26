C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      REAL FUNCTION FUNCMIN (PARM1, PARM2)
C=======================================================================

C   --*** FUNCMIN *** (ALGEBRA) Minimum function
C   --   Written by Amy Gilkey - revised 07/24/87
C   --
C   --FUNCMIN returns the minimum of two real numbers.  This is a function
C   --so that it can be passed to DOFNC2 (AMIN1 intrinsic cannot be passed).
C   --
C   --Parameters:
C   --   PARM1, PARM2 - IN - the minimum parameters

      FUNCMIN = MIN(PARM1, PARM2)

      RETURN
      END
