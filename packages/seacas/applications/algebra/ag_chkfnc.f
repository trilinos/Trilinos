C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHKFNC (FNAM, ENAM)
C=======================================================================

C   --*** CHKFNC *** (ALGEBRA) Check that the function names match
C   --   Written by Amy Gilkey - revised 07/24/87
C   --
C   --CHKFNC prints a program error message if the two function names
C   --do not match.  This is a debugging tool to check that the functions
C   --in EVAL are ordered correctly.
C   --
C   --Parameters:
C   --   FNAM - IN - the name of the function from the EVAL code
C   --   ENAM - IN - the name of the function from the equation entry

      CHARACTER*(*) FNAM, ENAM

      CHARACTER*80 ERRMSG

      IF (FNAM .NE. ENAM) THEN
         ERRMSG = 'Function ' // FNAM(:LENSTR(FNAM))
     &      // ' should be ' // ENAM(:LENSTR(ENAM))
         CALL PRTERR ('PROGRAM', ERRMSG(:LENSTR(ERRMSG)))
      END IF

      RETURN
      END
