C $Id: disctp.f,v 1.3 1991/03/22 16:07:18 gdsjaar Exp $
C $Log: disctp.f,v $
C Revision 1.3  1991/03/22 16:07:18  gdsjaar
C Fixed DATA statement screwed up in PI change
C
c Revision 1.2  1991/03/21  15:44:33  gdsjaar
c Changed all 3.14159... to atan2(0.0, -1.0)
c
c Revision 1.1.1.1  1990/11/30  11:06:09  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:07  gdsjaar
c Initial revision
c 
C
      LOGICAL FUNCTION DISCTP (ANGLE)
C***********************************************************************
C
C  FUNCTION DISCTP = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ANGLE IS
C                    WITHIN THE CURRENT DEFINITION OF A
C                    DISSECTION VERTEX
C
C***********************************************************************
C
      DATA EPS /.31/
C
      PI = ATAN2(0.0, -1.0)
      IF (ANGLE .GT. (PI + EPS)) THEN
         DISCTP=.TRUE.
      ELSE
         DISCTP=.FALSE.
      ENDIF
      RETURN
      END
