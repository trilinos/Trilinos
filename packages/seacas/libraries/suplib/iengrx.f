C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION IENGRX (IMAXE, IMINE)
C=======================================================================

C   --*** IENGX *** (STRLIB) Internal to NUMSTR
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --IENGRX returns the "best" engineering notation exponent for a
C   --minimum and maximum exponent.
C   --
C   --Parameters:
C   --   IMAXE - IN - the maximum exponent
C   --   IMINE - IN - the minimum exponent

      INTEGER IMAXE, IMINE

      IF (IMAXE .GT. 0) THEN
         IENGRX = INT ((IMAXE - 1) / 3) * 3
      ELSE
         IENGRX = INT ((IMAXE - 2) / 3) * 3
      END IF

      RETURN
      END
