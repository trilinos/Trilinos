C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: pack.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:06:37  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:29  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              pack
c
c DESCRIPTION:          Removes all spaces from an ASCII character
c                       string and counts the number of remaining
c                       non-space characters.
c
c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   subroutine
c
c USAGE:               call pack ( string, lstr )
c
c PARAMETERS:
c        character * (*) string    ( INPUT/OUTPUT )
c                       string to be stripped of spaces.
c        integer lstr               ( OUTPUT )
c                       number of non-space characters
c                       remaining.
c
c CALLS:                none
c
c GLOBAL VARIABLES REFERENCED:      none
c
c CALLING ROUTINE(S):      getins (BLOT)
c                          filhnd (BLOT)
c
c SYSTEM DEPENDENCIES:     none
c
c ======================================================================
c ======================================================================

      subroutine pack ( string, lstr )

c
c        parameters
c
      character * (*) string
      integer lstr

c
c        declarations
c
      integer length, pt1, pt2

c
c *************************************************************
c *************************************************************
c
c
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

c
c        fill remainder of string with spaces
c

      do 110 i = lstr + 1, length
         string(i:i) = ' '
  110 continue

      return
      end
