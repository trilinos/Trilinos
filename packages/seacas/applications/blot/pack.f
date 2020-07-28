C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================

c ROUTINE:              pack

c DESCRIPTION:          Removes all spaces from an ASCII character
c                       string and counts the number of remaining
c                       non-space characters.

c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511

c DATE:                 December 20, 1988

c TYPE OF SUBPROGRAM:   subroutine

c USAGE:               call pack ( string, lstr )

c PARAMETERS:
c        character * (*) string    ( INPUT/OUTPUT )
c                       string to be stripped of spaces.
c        integer lstr               ( OUTPUT )
c                       number of non-space characters
c                       remaining.

c CALLS:                none

c GLOBAL VARIABLES REFERENCED:      none

c CALLING ROUTINE(S):      getins (BLOT)
c                          filhnd (BLOT)

c SYSTEM DEPENDENCIES:     none

c ======================================================================
c ======================================================================

      subroutine pack ( string, lstr )

c        parameters

      character * (*) string
      integer lstr

c        declarations

      integer length, pt1, pt2

c *************************************************************
c *************************************************************

      length = len(string)

      pt1 = 0
      pt2 = 1
  100 continue
      if ( string(pt2:pt2) .ne. ' ' ) then
         pt1 = pt1 + 1
         string(pt1:pt1) = string(pt2:pt2)
      endif
      pt2 = pt2 + 1
      if ( pt2 .le. length ) go to 100
      lstr = pt1

c        fill remainder of string with spaces

      do 110 i = lstr + 1, length
         string(i:i) = ' '
  110 continue

      return
      end
