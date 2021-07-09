C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C************************************************************************
      subroutine wrprop(ndbi, ndbo, flag, nump, nval, names, values,
     *  namlen)
C************************************************************************
C     This subroutine with read properties from an input file and
C     write properties to an output text file
C     ndbi   - IN - file id to read from
C     ndbo   - IN - file id to write to
C     flag   - IN - flag indicating what type of properties to read/write
C     nump   - IN - number of properties
C     nval   - IN - number of values
C     names  - IN - names of the properties
C     values - IN - values of the properties

      integer ndbi, ndbo, flag, nump, nval
      character*(namlen) names(*)
      integer values(*)

      if (nump .eq. 0) return

      call exgpn (ndbi, flag, names, ierr)

      do 100 i = 1, nump
         call exgpa (ndbi, flag, names(i), values, ierr)
         write (ndbo,1000) '! Property Name: '
         write (ndbo,1000) names(i)
         write (ndbo,1000) '! Property Value(s): '
         write (ndbo,1010) (values(j), j = 1, nval)
 100  continue

 1000 format (A)
 1010 format (7(I10, 1X))
      return
      end
