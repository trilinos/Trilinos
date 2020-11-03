C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTRIM(LINE,L)
      CHARACTER*(*) LINE
      CHARACTER CH

      L1 = LEN(LINE)
 2240 IF (.NOT. (L1.GT.0)) GO TO 2260
      CH = LINE(L1:L1)
      IF (CH.EQ.' ') THEN
         GO TO 2250

      END IF

      GO TO 2260

 2250 L1 = L1 - 1
      GO TO 2240

 2260 CONTINUE
      L = L1
      RETURN

      END
