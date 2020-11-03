C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LOWFND (MXND, NUID, N, INDX, I, IOLD)
C***********************************************************************

C  SUBROUTINE LOWFND  =  LINEAR INDEXED SEARCH FOR MATCHING NUID VALUES

C***********************************************************************

      DIMENSION NUID(MXND), INDX(N)

      IBOT = 1
      ITOP = N

  100 CONTINUE
      II = (IBOT + ITOP)/2
      IF (NUID(INDX(II)) .EQ. NUID(I)) THEN
         IOLD = INDX(II)
         RETURN
      ELSE IF (NUID(INDX(II)) .GT. NUID(I)) THEN
         ITOP = II - 1
      ELSE
         IBOT = II + 1
      ENDIF
      IF (IBOT .LE. ITOP) GO TO 100

      IOLD = 0
      RETURN

      END
