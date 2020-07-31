C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c=======================================================================
      subroutine muntt (nelblk0, nelblk, nvarel, itt0, itt, istat)
C=======================================================================
C   --   ISTAT - IN - the status of each block:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n

      implicit none
      integer nelblk0, nelblk, nvarel
      integer istat(*)
      logical itt0 (nelblk0, nvarel)
      logical itt  (nelblk,  nvarel)

      integer ielb, ielb0, ivar

      ielb = 0
      do 10 ielb0 = 1, nelblk0
         if (istat(ielb0) .eq. 0) then
            ielb = ielb + 1
            do 20 ivar = 1, nvarel
               itt(ielb,ivar) = itt0(ielb0, ivar)
C               write (*,*) ielb0, ivar,itt0(ielb0, ivar)
 20         continue
         else if (istat(ielb0) .lt. 0) then
C ... Zero out the original truth table since we don't need to
C     read those variables off of the input database...
            do 30 ivar = 1, nvarel
               itt0(ielb0,ivar) = .false.
 30         continue

         else
            call prterr('PROGRAM',
     *       'Element block combination does not work if'
     *       //' there are element variables on the database.')
            stop
         end if
 10   continue
      if (ielb .ne. nelblk) then
         call prterr('PROGRAM', 'Problem in muntt')
         stop
      end if
      return
      end
