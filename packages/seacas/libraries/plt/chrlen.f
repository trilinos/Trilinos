C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION CHRLEN(S)
      CHARACTER*(*) S

      L = LEN(S)
      J = INDEX(S,CHAR(0))
      IF (J.EQ.0) THEN
         I = L

      ELSE

         I = J - 1
      END IF

   10 CONTINUE

      IF ((I.GT.0.AND.S(I:I).EQ.' ')) THEN
         I = I - 1
         GO TO 10

      END IF

      CHRLEN = I

      END
