C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SUBTRI (NPER, NEWPER, X, Y, NID, XSUB, YSUB, NIDSUB,
     &   M1, M2, IADD, ITRI, XCEN, YCEN)
C***********************************************************************

C  SUBROUTINE SUBPER = PUTS A SUBREGION'S PERIMETER INTO THE NPERIM
C                      ARRAYS

C***********************************************************************

      DIMENSION X (NPER), Y (NPER), NID (NPER)
      DIMENSION XSUB (NPER), YSUB (NPER), NIDSUB (NPER)

C  PUT SIDE ONE AND TWO INTO THE PERIMETER LIST

      KOUNT = 0
      DO 100 I = 1, M1 + M2 + 1
         KOUNT = KOUNT + 1
         J = I + IADD
         IF (J .GT. NPER)J = J - NPER
         XSUB (KOUNT) = X (J)
         YSUB (KOUNT) = Y (J)
         NIDSUB (KOUNT) = NID (J)
  100 CONTINUE

C  PUT SIDE THREE INTO THE LIST

      XDIF = XCEN - XSUB (KOUNT)
      YDIF = YCEN - YSUB (KOUNT)
      XINT = XDIF / DBLE(M1)
      YINT = YDIF / DBLE(M1)
      DO 110 I = 1, M1 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (ITRI * 100000) + M1 - I + 1
  110 CONTINUE

C  ENTER THE CENTER POINT

      KOUNT = KOUNT + 1
      XSUB (KOUNT) = XCEN
      YSUB (KOUNT) = YCEN
      NIDSUB (KOUNT) = 100000

C  PUT SIDE FOUR INTO THE LIST

      ITRI2 = ITRI + 2
      IF (ITRI2 .GT. 3)ITRI2 = ITRI2 - 3
      XDIF = X (IADD + 1) - XCEN
      YDIF = Y (IADD + 1) - YCEN
      XINT = XDIF / DBLE(M2)
      YINT = YDIF / DBLE(M2)
      DO 120 I = 1, M2 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (100000 * ITRI2) + I + 1
  120 CONTINUE

      NEWPER = KOUNT

      RETURN

      END
