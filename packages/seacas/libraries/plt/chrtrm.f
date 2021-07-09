C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHRTRM(LINE,L)
      CHARACTER*(*) LINE
      CHARACTER CH

      L1 = LEN(LINE)
 2180 IF (.NOT. (L1.GT.0)) GO TO 2200
      CH = LINE(L1:L1)
      IF (CH.EQ.' ') THEN
         GO TO 2190

      END IF

      IF (CH.EQ.CHAR(9)) THEN
         GO TO 2190

      END IF

      IF (CH.EQ.CHAR(0)) THEN
         GO TO 2190

      END IF

      GO TO 2200

 2190 L1 = L1 - 1
      GO TO 2180

 2200 CONTINUE
      L = L1
      RETURN

      END
