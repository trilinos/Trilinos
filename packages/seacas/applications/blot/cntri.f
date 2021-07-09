C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      REAL FUNCTION CNTRI (ICNTR)
C=======================================================================

C   --*** CNTRI *** (DETOUR) Return contour value
C   --   Written by Amy Gilkey - revised 06/04/87
C   --
C   --CNTRI returns the requested contour interval, which is either derived
C   --from CMIN and DELC or specified in CINTV.
C   --
C   --Parameters:
C   --   ICNTR - IN - the contour interval to be returned
C   --
C   --Common Variables:
C   --   Uses CINTOK, CMIN, DELC, CINTV of /CNTR/

      COMMON /CNTR/   CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC,
     &   CINTV(256), NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX
      LOGICAL CINTOK, LINCON, NOCMIN, NOCMAX

      IF (.NOT. CINTOK) THEN
         CNTRI = CMIN + (ICNTR-1) * DELC
      ELSE
         CNTRI = CINTV(ICNTR)
      END IF

      RETURN
      END
