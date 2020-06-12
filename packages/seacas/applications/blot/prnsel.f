C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: prnsel.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:07:54  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1992/05/22  22:37:29  gdsjaar
c Modified to handle more than 100,000 nodes/elements in printouts
c
c Revision 1.2  1990/12/14  08:55:18  gdsjaar
c Added RCS Id and Log to all files
c
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
