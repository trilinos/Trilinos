C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INDEXI_FQ (NMAX, IARRAY, N, INDX)
C***********************************************************************

C     INDEXI - SORT THE N ELEMENTS OF IARRAY WHOSE POSITION IS STORED
C              IN INDX.  ONLY THE ORDER OF THE INDEX ARRAY, INDX, IS
C              MODIFIED

C***********************************************************************

      DIMENSION IARRAY(NMAX), INDX(N)

      L = N/2 + 1
      IR = N
  100 CONTINUE
      IF (L .GT. 1) THEN
         L = L - 1
         INDXT = INDX(L)
         JQ = IARRAY(INDXT)
      ELSE
         INDXT = INDX(IR)
         JQ = IARRAY(INDXT)
         INDX(IR) = INDX(1)
         IR = IR - 1
         IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
            RETURN
         END IF
      END IF

      I = L
      J = L + L
  110 CONTINUE
      IF (J .LE. IR) THEN
         IF (J .LT. IR) THEN
            IF (IARRAY(INDX(J)) .LT. IARRAY(INDX(J + 1))) J = J + 1
         END IF
         IF (JQ .LT. IARRAY(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
         ELSE
            J = IR + 1
         END IF
         GO TO 110
      END IF
      INDX(I) = INDXT
      GO TO 100

      END
