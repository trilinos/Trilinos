C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      REAL FUNCTION MYSUM (PARM1, PARM2)
C=======================================================================

C   --*** MYSUM *** (ALGEBRA) Addition function
C   --   Written by Amy Gilkey - revised 07/24/87
C   --
C   --MYSUM returns the sum of two real numbers.  This is a function so that
C   --it can be passed to DOFNC2.
C   --
C   --Parameters:
C   --   PARM1, PARM2 - IN - the numbers to be summed

      MYSUM = PARM1 + PARM2

      RETURN
      END
