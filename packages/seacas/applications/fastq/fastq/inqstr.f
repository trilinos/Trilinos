C $Id: inqstr.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: inqstr.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1991/03/22 19:38:52  gdsjaar
C Fixed typo 0 was K0
C
c Revision 1.1.1.1  1990/11/30  11:10:02  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:10:00  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]UTIL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INQSTR (PROMPT, IANS)
C***********************************************************************
C
C  SUBROUTINE INQSTR = INPUTS CHARACTER STRINGS
C
C***********************************************************************
C
      CHARACTER* (*) PROMPT, IANS, HOLD*80
C
      IZ = 0
  100 CONTINUE
      CALL GETINP (IZ, IZ, PROMPT, HOLD, IOSTAT)
      IF (IOSTAT .EQ. 0) THEN
         CALL STRCUT (HOLD)
         IANS = HOLD (1:)
         RETURN
      ELSEIF (IOSTAT .LT. 0) THEN
         IANS = ' '
         RETURN
      ELSEIF (IOSTAT .GT. 0) THEN
         WRITE (*, 10010)
         GOTO 100
      ENDIF
C
10010 FORMAT (' BAD CHARACTER STRING  -  TRY AGAIN')
      END
