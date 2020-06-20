C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MATSTR (INSTR, MATCH, NLET)
C=======================================================================
C$Id: matstr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: matstr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:32  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:30  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:36  gdsjaar
c Initial revision
c

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
