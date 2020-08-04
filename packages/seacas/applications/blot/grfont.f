C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRFONT (IFNT)
C=======================================================================

C   --*** GRFONT *** (GRPLIB) Set font (PLT)
C   --   Written by Amy Gilkey - revised 01/22/87
C   --
C   --GRFONT sets a graphics font on the current device.
C   --
C   --Parameters:
C   --   IFNT - IN - the font type (as in /GRPCOM/)

C   --Routines Called:
C   --   PLTSTT - (PLTLIB) Set text parameter
C   --      7 = (KCHSPA) spacing between characters
C   --      12 = (KFONT) font (1=roman, 3=sanserif, 2=stick)

      PARAMETER (KCHSPA=7, KFONT=12)

      LOGICAL LDUM, PLTSTT

      IF (IFNT .EQ. 1) THEN
C      --stick font
         LDUM = PLTSTT (KCHSPA, 1.15)
         LDUM = PLTSTT (KFONT, 2.0)
      ELSE IF (IFNT .EQ. 2) THEN
C      --sanserif font
         LDUM = PLTSTT (KCHSPA, 1.15)
         LDUM = PLTSTT (KFONT, 3.0)
      ELSE IF (IFNT .EQ. 3) THEN
C      --roman font
         LDUM = PLTSTT (KCHSPA, 1.0)
         LDUM = PLTSTT (KFONT, 1.0)
      END IF

      RETURN
      END
