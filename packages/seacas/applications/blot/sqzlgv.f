C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SQZLGV (NPTIMS, IPTIMS, WHOTIM, VALIN, NPTOUT, VALOUT)
C=======================================================================

C   --*** SQZLGV *** (TPLOT) Compress time-dependent values
C   --   Written by Amy Gilkey - revised 11/06/87
C   --
C   --SQZLGV compresses time-dependent curve values so that only the
C   --values for whole time steps are in the curve.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   VALIN - IN - the input data
C   --   NPTOUT - OUT - the number of selected whole time steps
C   --   VALOUT - OUT - the output (compressed) data

      INTEGER IPTIMS(*)
      LOGICAL WHOTIM(*)
      REAL VALIN(*)
      REAL VALOUT(*)

      NPTOUT = 0
      DO 100 I = 1, NPTIMS
         IF (WHOTIM(IPTIMS(I))) THEN
            NPTOUT = NPTOUT + 1
            VALOUT(NPTOUT) = VALIN(I)
         END IF
  100 CONTINUE

      RETURN
      END
