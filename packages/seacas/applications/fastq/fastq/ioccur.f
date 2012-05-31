C $Id: ioccur.f,v 1.1 1990/11/30 11:10:28 gdsjaar Exp $
C $Log: ioccur.f,v $
C Revision 1.1  1990/11/30 11:10:28  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]IOCCUR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      FUNCTION IOCCUR (N, L, NEW)
C***********************************************************************
C
C  FUNCTION IOCCUR = CHECKS TO SEE IF NEW OCCURS IN  (L (I), I=1, N)
C
C***********************************************************************
C
C     RETURN 0 IF NEW DOES NOT OCCUR IN  (L (I), I=1, N)
C     RETURN 1 IF IT DOES
C
C***********************************************************************
C
      DIMENSION L (N)
C
      IF (N .LT. 1) THEN
         IOCCUR = 0
         RETURN
      ENDIF
C
      DO 100 I = 1, N
         IF (L (I) .EQ. NEW) THEN
            IOCCUR = 1
            RETURN
         ENDIF
  100 CONTINUE
      IOCCUR = 0
      RETURN
      END
