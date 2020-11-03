C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MINMXS (NPTSEL, IXSEL, PTS, VMIN, VMAX)
C=======================================================================

C   --*** MINMXS *** (ETCLIB) Calculate min/max value of selected points
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --MINMXS calculates the minimum and maximum of the selected data points.
C   --
C   --Parameters:
C   --   NPTSEL - IN - the number of selected points
C   --   IXSEL - IN - the node numbers of the selected nodes
C   --   PTS - IN - the points
C   --   VMIN, VMAX - OUT - the maximum and maximum value of the points

      INTEGER NPTSEL
      INTEGER IXSEL(*)
      REAL PTS(*)
      REAL VMIN, VMAX

      VMIN =  1.0E+30
      VMAX = -1.0E+30
      DO 10 IX = 1, NPTSEL
         I = IXSEL(IX)
         VMIN = MIN (VMIN, PTS(I))
         VMAX = MAX (VMAX, PTS(I))
   10 CONTINUE

      RETURN
      END
