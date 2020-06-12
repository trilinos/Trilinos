C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHOCMD (HEADER, LIST)
C=======================================================================
C$Id: shocmd.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: shocmd.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:16:19  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:16:18  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:45  gdsjaar
c Initial revision
c

C   --*** SHOCMD *** (ETCLIB) Display list of strings
C   --   Written by Amy Gilkey - revised 10/27/86
C   --
C   --SHOCMD displays the list of strings.  The heading for the list is
C   --dependent on HEADER:
C   --   COMMANDS - "Valid Commands"
C   --   null - no heading
C   --   other - HEADER
C   --
C   --Parameters:
C   --   HEADER - IN - the show option (see above)
C   --   LIST - IN - the string list, last entry must be ' '

C   --Routines Called:
C   --   LOCSTR - (STRLIB) Find string

      CHARACTER*(*) HEADER
      CHARACTER*(*) LIST(*)

      IF (HEADER .EQ. 'COMMANDS') THEN
         WRITE (*, 10) 'Valid Commands:'
      ELSE IF (HEADER .NE. ' ') THEN
         WRITE (*, 10) HEADER
      END IF

      NCMD = LOCSTR (' ', 999, LIST) - 1
      WRITE (*, 20) (LIST(I), I=1,NCMD)

      RETURN
   10 FORMAT (/, 1X, 5A)
   20 FORMAT ((4X, 7(A8, :, 2X)))
      END
