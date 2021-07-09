C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INBODY (MR, N9, IIN, IFOUND, IRPB, ADDOLD, NOROOM)
C***********************************************************************

C  SUBROUTINE INBODY = INPUTS A BODY LIST INTO THE DATABASE

C***********************************************************************

      DIMENSION IIN (IFOUND), IRPB (MR)

      LOGICAL NOROOM, ADDOLD

      NOROOM = .TRUE.
      IF (.NOT.ADDOLD)N9 = 0
      DO 120 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 130
         IF (N9 + 1 .GT. MR)RETURN

C  SEE IF THE REGION IS ALREADY IN THE BODY LIST

         DO 100 J = 1, N9
            IF (IRPB (J) .EQ. JJ)GOTO 110
  100    CONTINUE

         N9 = N9 + 1
         IRPB (N9) = JJ
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE

      NOROOM = .FALSE.
      RETURN

      END
