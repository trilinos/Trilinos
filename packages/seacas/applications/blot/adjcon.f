C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ADJCON (DELCOK)
C=======================================================================

C   --*** ADJCON *** (DETOUR) Adjust the contour parameters
C   --   Written by Amy Gilkey - revised 07/10/87
C   --
C   --ADJCON adjusts the contour parameters to be consistent with each
C   --other.  If CINTOK is false and the parameter DELCOK is true, the
C   --contour maximum (CMAX) is set; if DELCOK is false, the contour
C   --interval (DELC) is set.  If CINTOK is true, the contour minimum and
C   --maximum (CMIN and CMAX) and interval (DELC) are set.
C   --
C   --Parameters:
C   --   DELCOK - IN - true if CMAX is to be set, false if DELC is to be set
C   --
C   --Common Variables:
C   --   Uses CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC, CINTV of /CNTR/
C   --   Sets CMIN, CMAX, DELC of /CNTR/

      include 'cntr.blk'

      LOGICAL DELCOK

      IF (LINCON) THEN
         NC = NCNTR
      ELSE
         NC = NCNTR+1
      END IF

      IF (.NOT. CINTOK) THEN
         IF (DELCOK) THEN
            CMAX = CMIN + (NC-1) * DELC
         ELSE
            DELC = CMAX - CMIN
            IF (NC .GT. 1) DELC = (CMAX - CMIN) / (NC-1)
         END IF

      ELSE
         CMIN = CINTV(NC)
         CMAX = CINTV(NC)
         DELC = CMAX - CMIN
         IF (NC .GT. 1) DELC = (CMAX - CMIN) / (NC-1)
      END IF

      RETURN
      END
