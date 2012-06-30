C $Id: cornp.f,v 1.2 1991/03/21 15:44:30 gdsjaar Exp $
C $Log: cornp.f,v $
C Revision 1.2  1991/03/21 15:44:30  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:05:24  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:05:22  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]CORNP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION CORNP (ANGLE)
C***********************************************************************
C
C  FUNCTION CORNP = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ANGLE IS
C                   WITHIN THE CURRENT DEFINITION OF A CORNER
C
C***********************************************************************
C
      DATA EPS /.62/
C
      PI = ATAN2(0.0, -1.0)
      IF (ANGLE .LT. ( PI - EPS)) THEN
         CORNP=.TRUE.
      ELSE
         CORNP=.FALSE.
      ENDIF
C
      RETURN
C
      END
