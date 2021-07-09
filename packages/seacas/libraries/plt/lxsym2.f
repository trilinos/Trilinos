C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXSYM2(SYM,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) SYM,CH
      CHARACTER*65 S

      S =
     *'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_$.
     *'
      SYM = ' '
      CH = ILINE(JLINE:JLINE)
      KS = 0
 2560 IF (.NOT. (INDEX(S,CH).NE.0)) GO TO 2570
      CH = ILINE(JLINE:JLINE)
      JLINE = JLINE + 1
      KS = KS + 1
      SYM(KS:KS) = CH
      GO TO 2560

 2570 CONTINUE
      NS = KS
      LXSYM2 = KS .GT. 0
      RETURN

      END
