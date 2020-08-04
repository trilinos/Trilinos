C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INHOLE (MR, N7, N29, JPNTR, IIN, IFOUND, IFHOLE, NHPR,
     &   IHLIST, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE INHOLE  =  INPUTS A REGION'S HOLES INTO THE DATABASE

C***********************************************************************

      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION IIN(IFOUND)

      LOGICAL NOROOM, MERGE

      NOROOM = .TRUE.

C  ADD THE REGION INTO THE DATABASE

      J = JPNTR
      IFHOLE(J) = N29 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN(I)
         IF (JJ .EQ. 0) GO TO 110
         N29 = N29 + 1
         IF (N29 .GT. MR*2) RETURN
         IHLIST(N29) = JJ
  100 CONTINUE

  110 CONTINUE
      NHPR(J) = N29 - IFHOLE(J) + 1
      IF (NHPR(J) .LT. 1) THEN
         WRITE(*, 10000) J
         N29 = IFHOLE(J) - 1
      END IF

      NOROOM = .FALSE.
      RETURN

10000 FORMAT(' REGION:', I5, ' HAS LESS THAN ONE HOLE', /,
     &   ' THE HOLES FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')

      END
