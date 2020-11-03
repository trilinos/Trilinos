C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C************************************************************************
      subroutine rwpval1(ntxt, ndb, flag, nump, nval, id, names,
     &  values, *)
C************************************************************************
C     This subroutine with read properties from an input text file and
C     write properties to an output exodusII file
C     ntxt   - IN - file id to read from
C     ndb    - IN - file id to write to
C     flag   - IN - flag indicating what type of properties to read/write
C     nump   - IN - number of properties
C     nval   - IN - number of values
C     names  - IN - names of the properties
C     values - IN - values of the properties

      include 'exodusII.inc'

      integer ntxt, ndb, flag, nump, nval
      integer id(*)
      character*(mxstln) names(*)
      integer values(*)

      if (nump .eq. 0) return

      do 100 i = 1, nump
C ... Skip comment line
        read (ntxt, *, end=130, err=130)
        read (ntxt, '(A)', end=120, err=120) names(i)
C ... Skip comment line
        read (ntxt, *, end=130, err=130)
        read (ntxt, *, end=130, err=130) (values(j), j = 1, nval)

C ... Write the property (unless it is 'ID')
        if (names(i) .ne. 'ID                              ') then
          do 90 j = 1, nval
            call expp(ndb, flag, id(j), names(i), values(j), ierr)
 90       continue
        end if
 100  continue
      return

 120  continue
 130  continue
      CALL PRTERR('FATAL', 'Problem reading property')
      return 1
      end
