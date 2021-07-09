C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================

c ROUTINE:              Upcase

c DESCRIPTION:          Changes ASCII text strings to upper case.

c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511

c DATE:                 December 20, 1988

c TYPE OF SUBPROGRAM:   subroutine

c USAGE:               call upcase ( string )

c PARAMETERS:

c        character*(*) string -- ( Input and Output )
c                       String that is to be converted to
c                          upper case.

c CALLS:

c        len (INTRINSIC) -- returns length of character string.
c        ichar (INTRINSIC) -- returns ASCII integer value
c                             of passed character.
c        char (INTRINSIC) -- returns ASCII character assigned
c                            to passed integer value.

c GLOBAL VARIABLES REFERENCED:   none

c SYSTEM DEPENDENCIES:           none

c CALLING ROUTINE(S):            filhnd (BLOT)

c ======================================================================
c ======================================================================

      subroutine upcase_bl ( string )

c **********************************************************************

c        parameter

      character*(*) string

c **********************************************************************

c        declarations

      integer length
c          length of character string
      integer ccode
c          integer id of an ASCII character

c***********************************************************************
c***********************************************************************

      length = len(string)
      do 100 i = 1,length
         ccode = ichar(string(i:i))
         if ( ccode .ge. 97  .AND.  ccode .le. 122 )
     &      string(i:i) = char(ccode-32)
  100 continue

      return
      end
