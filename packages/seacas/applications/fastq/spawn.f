C $Id: spawn.f,v 1.2 1990/11/30 11:46:57 gdsjaar Exp $
C $Log: spawn.f,v $
C Revision 1.2  1990/11/30 11:46:57  gdsjaar
C Removed LIB$SPAWN call
C
c Revision 1.1.1.1  90/11/30  11:16:22  gdsjaar
c FASTQ Version 2.0X
c 
c Revision 1.1  90/11/30  11:16:21  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]SPAWN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SPAWN (VAXVMS)
C***********************************************************************
C
C  SUBROUTINE SPAWN = SPAWNS A PROCESS FOR ESCAPE OUT OF FASTQ
C
C***********************************************************************
C
C  VARIABLES USED:
C     VAXVMS = .TRUE. IF RUNNING ON A VAXVMS SYSTEM
C
C***********************************************************************
C
      LOGICAL VAXVMS
C
      IF (VAXVMS) THEN
         continue
      ELSE
         CALL MESAGE ('SPAWNING POSSIBLE ONLY ON VAXVMS SYSTEM')
      ENDIF
C
      END
