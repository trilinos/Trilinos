C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBERR (IOSTAT, ERRMSG)
C=======================================================================

C   --*** DBERR *** (EXOLIB) Display a database error message
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBERR displays a database read error message.
C   --
C   --Parameters:
C   --   IOSTAT - IN - the error status
C   --   ERRMSG - IN/OUT - the item being read when the error occurred;
C   --      compressed

C   --Routines Called:
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) ERRMSG

      CALL SQZSTR (ERRMSG, LSTR)
      WRITE (*, 10000) ERRMSG(:LSTR)
      IF (IOSTAT .LT. 0) THEN
         WRITE (*, 10010) IOSTAT, 'Unexpected end of file'
      ELSE IF (IOSTAT .EQ. 67) THEN
         WRITE (*, 10010) IOSTAT, 'Input record is too short'
      ELSE
         WRITE (*, 10010) IOSTAT
      END IF

      RETURN

10000  FORMAT (/, ' DATABASE ERROR - Reading ', A)
10010  FORMAT (3X, ' FORTRAN Error #', I6, :, ' - ', A)
      END
