C $Id: fexit.f,v 1.3 1999/01/27 15:17:47 gdsjaar Exp $
C $Log: fexit.f,v $
C Revision 1.3  1999/01/27 15:17:47  gdsjaar
C Added typical summary of mesh data on output.
C
C Better filename handling
C
C Cleaned up some character string handling
C
C Revision 1.2  1993/07/21 18:11:37  gdsjaar
C Removed message after stop
C
c Revision 1.1.1.1  1990/11/30  11:07:22  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:07:21  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]FEXIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &   TIME1, BATCH, VERSN)
C***********************************************************************
C
C  SUBROUTINE FEXIT = GRACEFULL FASTQ EXIT
C
C***********************************************************************
C
      CHARACTER*72 CIN (MCOM), VERSN*10, DATE*8, TIME*8
C
      LOGICAL IANS, WROTE, BATCH
C
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
C
      IF (.NOT.WROTE)THEN
         CALL MESAGE (' ')
         CALL MESAGE ('***********************************************')
         CALL MESAGE ('*  WARNING: NO OUTPUT FILE HAS BEEN WRITTEN   *')
         CALL MESAGE ('***********************************************')
         CALL INTRUP ('EXIT ANYWAY', IANS, MCOM, ICOM, JCOM, CIN, IIN,
     &      RIN, KIN)
         IF (.NOT.IANS)RETURN
      ENDIF
      CALL MESAGE (' ')
      CALL EXCPUS (TIME2)
      IF (BATCH)THEN
         CALL MESAGE ('FASTQ COMPLETED SUCCESSFULLY')
         CALL EXDATE (DATE)
         CALL EXTIME (TIME)
         WRITE (*, *)'             DATE: ', DATE
         WRITE (*, *)'             TIME: ', TIME
         WRITE (*, *)'          VERSION: ', VERSN
         WRITE (*, ' (A, I5)')'  CPU SECONDS USED: ', IFIX (TIME2-TIME1)
         CALL MESAGE (' ')
      ELSE
         WRITE (*, ' (A, I5)')' CPU SECONDS USED: ', IFIX (TIME2-TIME1)
         CALL MESAGE ('*--------------------------*')
         CALL MESAGE (' ')
         CALL VDESCP (10003, 0, 0)
         CALL PLTEND
      ENDIF
      stop
      END
