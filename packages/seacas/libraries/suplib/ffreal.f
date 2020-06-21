C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFREAL (IFLD, INTYP, RFIELD, EXPECT, DEFVAL, RVAL, *)
C=======================================================================
C$Id: ffreal.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffreal.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:38  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:37  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:27  gdsjaar
c Initial revision
c

C   --*** FFREAL *** (FFLIB) Parse free-field real
C   --   Written by Amy Gilkey - revised 02/24/86
C   --
C   --FFREAL parses a real field.  A default is supplied if the field
C   --is empty.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --   INTYP - IN - the input type from the free-field reader
C   --   RFIELD - IN - the real field
C   --   EXPECT - IN - the value to expect string, for error
C   --   DEFVAL - IN - the default value if field is empty
C   --   RVAL - OUT - the real value, set only if no error
C   --   * - return statement if the field is invalid; message is printed

      INTEGER IFLD
      INTEGER INTYP(*)
      REAL RFIELD(*)
      CHARACTER*(*) EXPECT
      REAL DEFVAL, RVAL

      CHARACTER*80 ERRMSG

      IF (INTYP(IFLD) .GE. 1) THEN
         RVAL = RFIELD(IFLD)
      ELSE IF (INTYP(IFLD) .LE. -1) THEN
         RVAL = DEFVAL
      ELSE
         ERRMSG = 'Expected ' // EXPECT
         CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         GOTO 100
      END IF

      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN

  100 CONTINUE
      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN 1
      END
