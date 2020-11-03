C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCSTR (STR, LENLST, STRLST)
C=======================================================================

C   --*** LOCSTR *** (STRLIB) Find string in list
C   --   Written by Amy Gilkey - revised 03/21/88
C   --
C   --LOCSTR returns the index of the given string in a list of strings.
C   --If the string is not in the list, returns 0.
C   --
C   --Parameters:
C   --   STR - IN - the string to be searched for
C   --   LENLST - IN - the number of strings in the list
C   --   STRLST - IN - the list of strings to be searched

      CHARACTER*(*) STR
      INTEGER LENLST
      CHARACTER*(*) STRLST(*)

      DO 10 LOCSTR = 1, LENLST
         IF (STR .EQ. STRLST(LOCSTR)) GOTO 20
   10 CONTINUE
      LOCSTR = 0

   20 CONTINUE
      RETURN
      END
