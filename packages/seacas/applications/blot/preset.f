C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRESET (LINPLT, QAPLT, DBORD0, DVIEW0)
C=======================================================================

C   --*** PRESET *** (BLOT) Initializes special graphics parameters
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --PRESET graphics parameters specific to the plot format:
C   --   o axes (set larger for line plots, same for mesh plots)
C   --   o tick marks (set larger for line plots, smaller for mesh plots)
C   --   o grid lines (set smaller for line plots)
C   --   o normal lines (set larger for line plots, same for mesh plots
C   --   o curve lines (set larger for line plots)
C   --
C   --It also sets the display area (dependent on QAPLT).
C   --
C   --Parameters:
C   --   LINPLT - IN - true if line versus mesh plot
C   --   QAPLT - IN - if true, set up as standard QA plot, else set up to
C   --      plot the entire screen (minus border)
C   --   DBORD0 - IN - the display and label area boundary
C   --   DVIEW0 - OUT - the display area boundary

C   --Routines Called:
C   --   PLTGTG - (PLTLIB) Get graph parameters (see PLTSTG)
C   --   PLTSTG - (PLTLIB) Set graph parameters:
C   --      28 = (KAXESZ) axis size
C   --      33 = (KTICSZ) tick-mark size
C   --      34 = (KGRISZ) grid-line size
C   --      29 = (KCRVSZ) line size
C   --      15 = (KONGRI) no grid
C   --      37 = (KZLINE) delete dashed line at 0

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (KAXESZ=28, KTICSZ=33, KGRISZ=34, KCRVSZ=29)
      PARAMETER (KZLINE=37)
      PARAMETER (KONGRI=15)

      LOGICAL LINPLT, QAPLT
      REAL DBORD0(KTOP), DVIEW0(KTOP)

      LOGICAL LDUM, PLTGTG, PLTSTG

      LOGICAL FIRST
      SAVE FIRST
      SAVE SZAXES, SZTICK, SZGRID, SZCURV

      DATA FIRST / .TRUE. /

      IF (QAPLT) THEN
         DVIEW0(KLFT) = DBORD0(KLFT) + 0.09
         DVIEW0(KRGT) = DVIEW0(KLFT) + 0.62
         DVIEW0(KBOT) = DBORD0(KBOT) + 0.10
         DVIEW0(KTOP) = DVIEW0(KBOT) + 0.62
      ELSE
         BORDER = 0.1
         DVIEW0(KLFT) = DBORD0(KLFT) + BORDER
         DVIEW0(KRGT) = DBORD0(KRGT) - BORDER
         DVIEW0(KBOT) = DBORD0(KBOT) + BORDER
         DVIEW0(KTOP) = DBORD0(KTOP) - BORDER
      END IF

C   --Save default device settings

      IF (FIRST) THEN
         LDUM = PLTGTG (KAXESZ, SZAXES)
         LDUM = PLTGTG (KTICSZ, SZTICK)
         LDUM = PLTGTG (KGRISZ, SZGRID)
         LDUM = PLTGTG (KCRVSZ, SZCURV)
         FIRST = .FALSE.
      END IF

      IF (LINPLT) THEN
         LDUM = PLTSTG (KAXESZ, 1.75*SZAXES)
         LDUM = PLTSTG (KTICSZ, 1.75*SZTICK)
         LDUM = PLTSTG (KGRISZ, 0.50*SZGRID)
         LDUM = PLTSTG (KCRVSZ, 1.50*SZCURV)
      ELSE
         LDUM = PLTSTG (KAXESZ, SZAXES)
         LDUM = PLTSTG (KTICSZ, 0.50*SZTICK)
      END IF

C   --Set no line at zero
      LDUM = PLTSTG (KZLINE, 0.0)

C   --Set no grid
      LDUM = PLTSTG (KONGRI, 0.0)

      RETURN
      END
