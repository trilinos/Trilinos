C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFCHAR (IFLD, INTYP, CFIELD, DEFVAL, CVAL)
C=======================================================================
C$Id: ffchar.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffchar.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:20  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:19  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:24  gdsjaar
c Initial revision
c

C   --*** FFCHAR *** (FFLIB) Parse free-field character string
C   --   Written by Amy Gilkey - revised 02/24/86
C   --
C   --FFCHAR parses a character field.  A default is supplied if the
C   --field is empty.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --   INTYP - IN - the input type from the free-field reader
C   --   CFIELD - IN - the character fields
C   --   DEFVAL - IN - the default value if field is empty
C   --   CVAL - OUT - the character value

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) DEFVAL, CVAL

      IF (INTYP(IFLD) .GE. 0) THEN
         CVAL = CFIELD(IFLD)
      ELSE IF (INTYP(IFLD) .LE. -1) THEN
         CVAL = DEFVAL
      END IF

      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN
      END
