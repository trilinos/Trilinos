C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      subroutine screrr(iunit, filin, lnam, name1, name2)
      integer iunit, lnam
      character*(*) filin, name1, name2
      character*32 tmpstr

      write(tmpstr, *) 'fort.',iunit
      call pckstr(1,tmpstr)

      if (filin(:lnam) .eq. tmpstr(:lnam)) then
         write (*,100) tmpstr(:lnam), name1, name2
 100     format(
     $'***********************************************************',/,
     $'***********************************************************',/,
     $' NOTE: The specified database is "',A,'"',/
     $'       This is usually the result of incorrectly running',/,
     $'       $ACCESS/bin/',A,' instead of $ACCESS/etc/',A,'.',//,
     $'       Please check your path and make sure that',/,
     $'       $ACCESS/etc precedes $ACCESS/bin.',/,
     $'***********************************************************',/,
     $'***********************************************************')
      end if
      return
      end
