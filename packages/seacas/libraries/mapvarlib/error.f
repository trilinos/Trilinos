C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,ERROR
      SUBROUTINE ERROR (SUBNAM,MESSAGE,LABEL1,I,LABEL2,J,LABEL3,WORD,
     1  ISTOP)

C     ******************************************************************

C     SUBROUTINE TO PRINT ERROR MESSAGE AND TERMINATE EXECUTION

C     Calls subroutine CLSFIL

C     Called by everything

C     ******************************************************************

      CHARACTER*(*) SUBNAM,MESSAGE,LABEL1,LABEL2,LABEL3,WORD

      include 'tapes.blk'

C     ******************************************************************

      WRITE (NOUT, 60)
      WRITE (NTPOUT, 60)
      WRITE (NOUT, 10) SUBNAM
      WRITE (NTPOUT, 10) SUBNAM
      WRITE (NOUT, 20) MESSAGE
      WRITE (NTPOUT, 20) MESSAGE
      WRITE (NOUT, 30)
      WRITE (NTPOUT, 30)
      IF (LABEL1.NE.' ') THEN
        WRITE (NOUT, 40) LABEL1,I
        WRITE (NTPOUT, 40) LABEL1,I
      END IF
      IF (LABEL2.NE.' ') THEN
        WRITE (NOUT, 40) LABEL2,J
        WRITE (NTPOUT, 40) LABEL2,J
      END IF
      IF (LABEL3.NE.' ') THEN
        WRITE (NOUT, 50) LABEL3,WORD
        WRITE (NTPOUT, 50) LABEL3,WORD
      END IF
      WRITE (NOUT, 60)
      WRITE (NTPOUT, 60)

      IF (ISTOP.EQ.0) RETURN

      CALL CLSFIL

      STOP 'ERROR'

   10 FORMAT (/,10X,' ERROR FOUND IN - ' ,A)
   20 FORMAT (/,10X,' DESCRIPTION - ' ,A)
   30 FORMAT (/,10X,' RELEVANT PARAMETERS - ')
   40 FORMAT (/,15X, A, ' = ' ,I10)
   50 FORMAT (/,15X, A, ' = ', A)
   60 FORMAT (/,10X,'* * * * * * * * * * * * * * * * * * * * * * * * *',
     *' * * * * * ',/,10X,'* * * * * * * * * * * * * * * * * * * * * *',
     *' * * * * * * * *',/)
      END
