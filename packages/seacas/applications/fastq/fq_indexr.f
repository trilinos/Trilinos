C $Id: indexr.f,v 1.1 1990/11/30 11:09:29 gdsjaar Exp $
C $Log: indexr.f,v $
C Revision 1.1  1990/11/30 11:09:29  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]INDEXR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INDEXR (NMAX, RARRAY, N, INDX)
C***********************************************************************
C
C     INDEXI - SORT THE N ELEMENTS OF RARRAY WHOSE POSITION IS STORED
C              IN INDX.  ONLY THE ORDER OF THE INDEX ARRAY, INDX, IS
C              MODIFIED
C
C***********************************************************************
C
      DIMENSION RARRAY(NMAX), INDX(N)
C
      L = N/2 + 1
      IR = N
  100 CONTINUE
      IF (L .GT. 1) THEN
         L = L - 1
         INDXT = INDX(L)
         Q = RARRAY(INDXT)
      ELSE
         INDXT = INDX(IR)
         Q = RARRAY(INDXT)
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
            IF (RARRAY(INDX(J)) .LT. RARRAY(INDX(J + 1))) J = J + 1
         END IF
         IF (Q .LT. RARRAY(INDX(J))) THEN
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
