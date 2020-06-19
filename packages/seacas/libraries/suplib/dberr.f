C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBERR (IOSTAT, ERRMSG)
C=======================================================================
C$Id: dberr.f,v 1.4 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dberr.f,v $
CRevision 1.4  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.3  1991/09/30 20:08:33  gdsjaar
CIncreased error number format for Cray
C
c Revision 1.2  1991/02/04  08:34:36  gdsjaar
c Changed IOSTAT format to I3
c
c Revision 1.1.1.1  90/08/14  16:12:31  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:12:30  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:08  gdsjaar
c Initial revision
c

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
