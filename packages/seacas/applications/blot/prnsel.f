C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRNSEL (NSEL, NTOT, VALNAM)
C=======================================================================

C   --*** PRNSEL *** (BLOT) Print number of values selected
C   --   Written by Amy Gilkey - revised 03/02/88
C   --
C   --PRNSEL prints the number of values selected.
C   --
C   --Parameters:
C   --   NSEL - IN - the number of values selected
C   --   NTOT - IN - the total number of values
C   --   VALNAM - IN - the name of the value being checked (plural)

      CHARACTER*(*) VALNAM

      CHARACTER*80 STRING

      IF (NSEL .LE. 0) THEN
         WRITE (STRING, '(5A)') 'No ', VALNAM, ' are selected'
         CALL PRTERR ('CMDWARN', STRING(:LENSTR(STRING)))
      ELSE
         WRITE (STRING, 10000) NSEL, NTOT, VALNAM
10000     FORMAT (I6, ' of ', I6, ' ', A, ' selected')
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
      END IF

      RETURN
10010  FORMAT (4X, A)
      END
