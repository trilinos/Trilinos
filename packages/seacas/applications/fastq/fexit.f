C    Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &   TIME1, BATCH, VERSN)
C***********************************************************************

C  SUBROUTINE FEXIT = GRACEFUL FASTQ EXIT

C***********************************************************************

      CHARACTER*72 CIN (MCOM), VERSN*10, DATE*8, TIME*8

      LOGICAL IANS, WROTE, BATCH

      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)

      call addlog (VERSN(:LENSTR(VERSN)))
      IF (.NOT.WROTE)THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('***********************************************')
         CALL MESSAGE('*  WARNING: NO OUTPUT FILE HAS BEEN WRITTEN   *')
         CALL MESSAGE('***********************************************')
         CALL INTRUP ('EXIT ANYWAY', IANS, MCOM, ICOM, JCOM, CIN, IIN,
     &      RIN, KIN)
         IF (.NOT.IANS)RETURN
      ENDIF
      CALL MESSAGE(' ')
      CALL EXCPUS (TIME2)
      IF (BATCH)THEN
         CALL MESSAGE('FASTQ COMPLETED SUCCESSFULLY')
         CALL EXDATE (DATE)
         CALL EXTIME (TIME)
         WRITE (*, *)'             DATE: ', DATE
         WRITE (*, *)'             TIME: ', TIME
         WRITE (*, *)'          VERSION: ', VERSN
         WRITE (*, ' (A, I5)')'  CPU SECONDS USED: ', INT (TIME2-TIME1)
         CALL MESSAGE(' ')
      ELSE
         WRITE (*, ' (A, I5)')' CPU SECONDS USED: ', INT (TIME2-TIME1)
         CALL MESSAGE('*--------------------------*')
         CALL MESSAGE(' ')
         CALL VDESCP (10003, 0, 0)
         CALL PLTEND
      ENDIF
      call mdfree()
      stop
      END
