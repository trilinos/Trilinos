C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE SELBLK (NUMSEL, IXSEL, NELBLK, LISBLK, NUMELB, NUMLNK,
     *  LINK, ISCR, NUMNP, EBTYPE)
C=======================================================================

      include 'exodusII.inc'

      INTEGER IXSEL(*)
      INTEGER LISBLK(0:*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      INTEGER ISCR(*)
      LOGICAL SELECT
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*40 STRA
      CHARACTER*132 MSG

      do 80 i=1, numnp
        iscr(i) = 0
 80   continue

      islnk = 1
      do 100 ielb = 1, nelblk
        select = .false.
        do 90 ix = 1, lisblk(0)
          if (ielb .eq. lisblk(ix)) then
            select = .true.
          end if
 90     continue
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
          numnod = numlnk(ielb)
        else
          numnod = numlnk(ielb) * numelb(ielb)
        end if
        if (select) then
          call selblk1(ielb, numnod, link(islnk), iscr)
        end if
        ISLNK = ISLNK + numnod
 100  CONTINUE

      numsel = 0
      do 120 i=1, numnp
        if (iscr(i) .gt. 0) then
          numsel = numsel + 1
          ixsel(numsel) = i
        end if
 120  continue

      write (stra, 10000) numsel
      call pckstr(1, stra)
      MSG = STRA(:lenstr(stra)) // ' nodes selected'
      call prterr('CMDSPEC', MSG)
10000 format(I12)
      return
      end

      subroutine selblk1(ielb, numnod, link, iscr)

      integer link(*)
      integer iscr(*)

      do i=1, numnod
        node = link(i)
        iscr(node) = iscr(node) + 1
      end do
      return
      end
