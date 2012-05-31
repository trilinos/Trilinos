C $Id: indexi.f,v 1.1 1990/11/30 11:09:26 gdsjaar Exp $
C $Log: indexi.f,v $
C Revision 1.1  1990/11/30 11:09:26  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]INDEXI.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INDEXI (NMAX, IARRAY, N, INDX)
C***********************************************************************
C
C     INDEXI - SORT THE N ELEMENTS OF IARRAY WHOSE POSITION IS STORED
C              IN INDX.  ONLY THE ORDER OF THE INDEX ARRAY, INDX, IS
C              MODIFIED
C
C***********************************************************************
C
      DIMENSION IARRAY(NMAX), INDX(N)
C
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
C
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
C
      END
