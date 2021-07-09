C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MATSTR (INSTR, MATCH, NLET)
C=======================================================================

C   --*** MATSTR *** (STRLIB) Check if string matches
C   --   Written by Amy Gilkey - revised 07/01/87
C   --
C   --MATSTR true iff the input string is equal to the match string.
C   --Only NLET letters must be in the input string to match, but if more
C   --letters are given, they must match the match string exactly.
C   --
C   --Parameters:
C   --   INSTR - IN - the input string
C   --   MATCH - IN - the match string
C   --   NLET - IN - number of letters that must match; 0 for exact match

      CHARACTER*(*) INSTR
      CHARACTER*(*) MATCH
      INTEGER NLET

      IF (NLET .LE. 0) THEN
         MATSTR = INSTR .EQ. MATCH
      ELSE
         LMATCH = LENSTR (MATCH)
         LMIN = MIN (LMATCH, NLET)
         LINSTR = LENSTR (INSTR)
         IF ((LINSTR .LE. LMATCH) .AND. (LINSTR .GE. LMIN)) THEN
            IF (LMIN .LT. LINSTR) LMIN = LINSTR
            MATSTR = INSTR(:LMIN) .EQ. MATCH(:LMIN)
         ELSE
            MATSTR = .FALSE.
         END IF
      END IF

      RETURN
      END
