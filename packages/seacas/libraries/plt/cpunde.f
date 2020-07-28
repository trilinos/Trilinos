C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPUNDE(IU)
      CHARACTER*40 CPUHLB
      COMMON /CPHLBN/CPUHLB
      INTEGER CPUNIT(20)
      INTEGER CPUNIF
      COMMON /CPUN/CPUNIT,CPUNIF
      CHARACTER*20 TERM
      CHARACTER*80 LINBUF(24)
      COMMON /REBUF/LINBUF,TERM
      INTEGER LINLEN(24)
      COMMON /REBUF2/NL,LINLEN

      J = IU - 80 + 1
      CPUNIT(J) = IU
      RETURN

      END
