C $Id: lowfnd.f,v 1.1 1990/11/30 11:11:28 gdsjaar Exp $
C $Log: lowfnd.f,v $
C Revision 1.1  1990/11/30 11:11:28  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]LOWFND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LOWFND (MXND, NUID, N, INDX, I, IOLD)
C***********************************************************************
C
C  SUBROUTINE LOWFND  =  LINEAR INDEXED SEARCH FOR MATCHING NUID VALUES
C
C***********************************************************************
C
      DIMENSION NUID(MXND), INDX(N)
C
      IBOT = 1
      ITOP = N
C
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
C
      IOLD = 0
      RETURN
C
      END
