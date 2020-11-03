C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRPBEG
C=======================================================================

C   --*** GRPBEG *** (GRPLIB) Begin a plot by setting defaults (PLT)
C   --   Written by Amy Gilkey - revised 03/01/88
C   --
C   --GRPBEG begins a plot by enabling the interrupt flag, clearing the
C   --screen, and setting standard graphics parameters (such as character
C   --sizes).  The graphics parameters are saved the first time through
C   --this routine and reset to the saved values on following calls.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Initialize interrupt
C   --   PLTBGN - (PLTLIB) Erase display surface
C   --   PLTGTG - (PLTLIB) Get graph parameter (see PLTSTG)
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      22, 47 = (KXNUMS, KYNUMS) X, Y axis number size
C   --      23, 48 = (KXLABS, KYLABS) X, Y axis label size
C   --   GRCOLR - (GRPLIB) Set color
C   --   GRSNAP - (GRPLIB) Handle device frame snapping

      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)
      PARAMETER (KTICSZ=33)

      LOGICAL CPUIFC
      LOGICAL LDUM, PLTGTG, PLTSTG

      LOGICAL FIRST
      SAVE FIRST
      SAVE SZNUM, SZLAB, SZTIC

      DATA FIRST / .TRUE. /

C   --Enable interrupt flag
      LDUM = CPUIFC (.TRUE.)

C   --Erase display surface
      CALL PLTBGN

C   --Start plot for snapping
      CALL GRSNAP ('START', 0)

C   --Save standard graphics parameters
      IF (FIRST) THEN
         LDUM = PLTGTG (KXNUMS, SZNUM)
         LDUM = PLTGTG (KXLABS, SZLAB)
         LDUM = PLTGTG (KTICSZ, SZTIC)
         FIRST = .FALSE.
      END IF

C   --Set up label/numbering size
      LDUM = PLTSTG (KXNUMS, SZNUM)
      LDUM = PLTSTG (KXLABS, SZLAB)
      LDUM = PLTSTG (KYNUMS, SZNUM)
      LDUM = PLTSTG (KYLABS, SZLAB)
      LDUM = PLTSTG (KTICSZ, SZTIC)

C   --Set to standard color
      CALL GRCOLR (0)

      RETURN
      END
