C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MEMTRC(MEMRY)
      INTEGER MEMRY(*)

      IKL = MEMRY(1)
      IKM = MEMRY(2)
      WRITE (6,*) IKL,IKM
      IP = 3
 2800 CONTINUE
      IPN = MEMRY(ABS(IP))
      PRINT *,IP,IPN
      IP = IPN
      IF (.NOT. (IPN.EQ.0)) GO TO 2800
      RETURN

      END
