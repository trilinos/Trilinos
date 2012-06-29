C $Id: inhole.f,v 1.1 1990/11/30 11:09:43 gdsjaar Exp $
C $Log: inhole.f,v $
C Revision 1.1  1990/11/30 11:09:43  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INHOLE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INHOLE (MR, N7, N29, JPNTR, IIN, IFOUND, IFHOLE, NHPR,
     &   IHLIST, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INHOLE  =  INPUTS A REGION'S HOLES INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION IIN(IFOUND)
C
      LOGICAL NOROOM, MERGE
C
      NOROOM = .TRUE.
C
C  ADD THE REGION INTO THE DATABASE
C
      J = JPNTR
      IFHOLE(J) = N29 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN(I)
         IF (JJ .EQ. 0) GO TO 110
         N29 = N29 + 1
         IF (N29 .GT. MR*2) RETURN
         IHLIST(N29) = JJ
  100 CONTINUE
C
  110 CONTINUE
      NHPR(J) = N29 - IFHOLE(J) + 1
      IF (NHPR(J) .LT. 1) THEN
         WRITE(*, 10000) J
         N29 = IFHOLE(J) - 1
      END IF
C
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT(' REGION:', I5, ' HAS LESS THAN ONE HOLE', /,
     &   ' THE HOLES FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')
C
      END
