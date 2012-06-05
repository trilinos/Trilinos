C $Id: inbody.f,v 1.1 1990/11/30 11:09:17 gdsjaar Exp $
C $Log: inbody.f,v $
C Revision 1.1  1990/11/30 11:09:17  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INBODY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INBODY (MR, N9, IIN, IFOUND, IRPB, ADDOLD, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INBODY = INPUTS A BODY LIST INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IIN (IFOUND), IRPB (MR)
C
      LOGICAL NOROOM, ADDOLD
C
      NOROOM = .TRUE.
      IF (.NOT.ADDOLD)N9 = 0
      DO 120 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 130
         IF (N9 + 1 .GT. MR)RETURN
C
C  SEE IF THE REGION IS ALREADY IN THE BODY LIST
C
         DO 100 J = 1, N9
            IF (IRPB (J) .EQ. JJ)GOTO 110
  100    CONTINUE
C
         N9 = N9 + 1
         IRPB (N9) = JJ
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
C
      NOROOM = .FALSE.
      RETURN
C
      END
