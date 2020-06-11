C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: upcase.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:17:19  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:59:25  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              Upcase
c
c DESCRIPTION:          Changes ASCII text strings to upper case.
c
c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   subroutine
c
c USAGE:               call upcase ( string )
c
c PARAMETERS:
c
c        character*(*) string -- ( Input and Output )
c                       String that is to be converted to
c                          upper case.
c
c CALLS:
c
c        len (INTRINSIC) -- returns length of character string.
c        ichar (INTRINSIC) -- returns ASCII integer value
c                             of passed character.
c        char (INTRINSIC) -- returns ASCII character assigned
c                            to passed integer value.
c
c
c GLOBAL VARIABLES REFERENCED:   none
c
c SYSTEM DEPENDENCIES:           none
c
c CALLING ROUTINE(S):            filhnd (BLOT)
c
c ======================================================================
c ======================================================================
c
      subroutine upcase_bl ( string )
c
c **********************************************************************
c
c        parameter
c
      character*(*) string
c
c **********************************************************************
c
c        declarations
c
      integer length
c          length of character string
      integer ccode
c          integer id of an ASCII character
c
c***********************************************************************
c***********************************************************************
c
      length = len(string)
      do 100 i = 1,length
         ccode = ichar(string(i:i))
         if ( ccode .ge. 97  .AND.  ccode .le. 122 )
     &      string(i:i) = char(ccode-32)
  100 continue

      return
      end
