C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHRRVC(BUFF,TXT,L)
      CHARACTER*16 LOCTXT
      CHARACTER*(*) TXT

      LT = LEN(TXT)
      WRITE (LOCTXT,'(1pg13.6)') BUFF
      DO 2110 I = 1,LT
         IF (LOCTXT(I:I).NE.' ') THEN
            L1 = I
            GO TO 2120

         END IF

 2110 CONTINUE
 2120 CONTINUE
      DO 2130 I = L1,LT
         IF (LOCTXT(I:I).EQ.' ') THEN
            L2 = I
            GO TO 2140

         END IF

 2130 CONTINUE
 2140 CONTINUE
      TXT = LOCTXT(L1:L2)
      L = L2 - L1
      RETURN

      END
