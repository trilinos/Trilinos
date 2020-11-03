C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EXPMAX (LABSID, VMIN, VMAX)
C=======================================================================

C   --*** EXPMAX *** (BLOT) Expands min/max values
C   --   Written by Amy Gilkey - revised 10/07/87
C   --
C   --EXPMAX expands the minimum and maximum values by 2.5% each.
C   --It also expands the appropriate limit by 2.5% for numbering.
C   --
C   --Parameters:
C   --   LABSID - IN - if 'FIRST', expand for numbering on left side;
C   --      if 'LAST', expand for numbering on right side
C   --   VMIN/VMAX - IN/OUT - the minimum and maximum axis values

      CHARACTER*(*) LABSID
      REAL VMIN, VMAX

      IF (VMIN .EQ. VMAX) THEN
         IF (VMIN .EQ. 0.0) THEN
            VMIN = -1.0
            VMAX = 1.0
         ELSE
            VRNG = 0.025 * ABS(VMIN)
            VMIN = VMIN - VRNG
            VMAX = VMAX + VRNG
         END IF

      ELSE
         VRNG = 0.025 * (VMAX - VMIN)
         VMIN = VMIN - VRNG
         VMAX = VMAX + VRNG
         IF (LABSID .EQ. 'FIRST') VMIN = VMIN - 0.025 * (VMAX - VMIN)
         IF (LABSID .EQ. 'LAST') VMAX = VMAX + 0.025 * (VMAX - VMIN)
      END IF

      RETURN
      END
