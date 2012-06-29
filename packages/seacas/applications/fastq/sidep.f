C $Id: sidep.f,v 1.2 1991/03/21 15:45:18 gdsjaar Exp $
C $Log: sidep.f,v $
C Revision 1.2  1991/03/21 15:45:18  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:15:43  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:15:41  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]SIDEP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION SIDEP (ANGLE)
C***********************************************************************
C
C  FUNCTION SIDEP = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ANGLE IS
C                   WITHIN THE CURRENT DEFINITION OF A SIDE
C
C***********************************************************************
C
      DATA EPS /1.27/
C
      PI = ATAN2(0.0, -1.0)
      IF ( (ANGLE .GT. (PI - EPS)) .AND. (ANGLE .LT. (PI + EPS)) ) THEN
         SIDEP=.TRUE.
      ELSE
         SIDEP=.FALSE.
      ENDIF
      RETURN
C
      END
